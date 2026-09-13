// File: ChargeDensity/Internal/FieldMixer.C  The G-SPACE mixing vocabulary shared by Kerker, Pulay and the
// polarized composition: the mixable FIELD (GField), the field-mixer and history faces, the joint
// extrapolation protocol, the raster-space Kerker step and the ΔG_Map arithmetic.  Internal: the only
// consumers are the concrete mixers (and their unit tests).
module;
#include <memory>
#include <complex>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <deque>
#include <vector>
export module qchem.ChargeDensity.Internal.FieldMixer;
export import qchem.ChargeDensity.FourierDensity;         // FourierDensity, ΔG_Map
export import qchem.ChargeDensity.FourierMixCD;           // FourierMixCD (the field's PRESENTATION as a density)
import qchem.BasisSet.G_FieldEvaluator;            // G_SpectralFilter: the raster Kerker step (0.5(f2))
import qchem.Math.DIIS;                             // the shared Pulay/DIIS bordered-solve engine
import qchem.Blaze;                                 // rsmat_t/rvec_t/ivec3_t + blazem::zero
import qchem.Types;

export namespace qchem::ChargeDensity
{

// --- The RASTER-SPACE Kerker step (the raw-XC shadow of FourierMixCD::KerkerMix; doc/GPWPlan 0.5(f2)) ---
// rho_mix(r) = rho_in + alpha * F^-1[ G^2/(G^2+G0^2) F(rho_out - rho_in) ] over the FULL box.  The filter is
// a SMOOTH multiplier (no truncation -> no Gibbs), k(0)=0 conserves charge -- identical algebra to the ball
// mix but over every mode the raster represents, so the raw feed's out-of-ball content keeps mixing at k~1.
inline rvec_t RasterKerker(const BasisSet::G_SpectralFilter& ge, const rvec_t& in, const rvec_t& out,
                           double alpha, double G0)
{
    const double G0sq=G0*G0;
    rvec_t delta=out; delta-=in;
    rvec_t mix=in;
    mix+=alpha*ge.ApplySpectralFilter(delta, [G0sq](double g2){return g2/(g2+G0sq);});
    return mix;
}

//! THE MIXABLE ENTITY: a G-space field \f$\tilde f(G)\f$ plus its optional raw-raster shadow.
//!
//! This type exists to make \f$m=\rho_\uparrow-\rho_\downarrow\f$ a FIRST-CLASS SUBJECT OF MIXING (user,
//! 2026-08-08).  A ρ̃ mixer's real precondition was never "a density" -- Kerker and Pulay only ever touched
//! the \c FourierDensity face -- yet they DECLARED \c Mix(cd_t&,...), i.e. a density-matrix-backed density.
//! That wider-than-necessary precondition is precisely what made m unsubstitutable: m has no density matrix,
//! no charge and no positivity, so it could never satisfy a requirement none of those mixers actually had.
//! Naming the true subject fixes it -- ρ and m are now THE SAME TYPE, and any field mixer takes either.
//! (Faking a density instead -- a \c tDM_CD with asserting stubs -- is the anti-pattern this codebase
//! avoids: capabilities live only on the types that have them.)
//!
//! NB the asymmetry is real and deliberate: the mix INPUT is a field (ρ or m equally), while the Fock
//! OUTPUT is genuinely a density -- the Hamiltonian consumes ρ, never m.  So substitutability belongs on
//! the input side, which is where it now lives.
struct GField
{
    ΔG_Map tilde;    //!< \f$\tilde f(G)\f$ on the fit basis's difference set
    rvec_t raster;   //!< the raw-raster shadow (doc/GPWPlan 0.5(f2)); empty = the raw pipeline is off
};

class tFieldExtrapolator;   // the HISTORY face, below

//! A mixer whose subject is a \c GField -- i.e. one that works purely in G space (Kerker, Pulay; NOT the
//! DM-based linear mixer, which genuinely needs the density matrix).  \c tDensityMixer::Mix is then just the
//! adapter that extracts the field from a density and calls this.
class tFieldMixer
{
public:
    virtual ~tFieldMixer() {}
    //! Fold the freshly collocated \a out into the running mixed field; returns the residual ‖out−in‖_∞.
    virtual double MixField(const GField& out) = 0;
    //! The running mixed density -- the PRESENTATION of the mixed field (what \c FockDensity returns, and
    //! what a caller recombining spin channels reads ρ̃ back from).
    virtual const FourierMixCD& Mixed() const = 0;
    //! \brief DOES THIS MIXER CARRY HISTORY? -- the property that decides whether a multi-channel caller may
    //! run one of these PER CHANNEL.  Returns the staging face (below), or \c nullptr for a memoryless filter.
    //!
    //! A memoryless FILTER (Kerker's \f$G^2/(G^2+G_0^2)\f$) is linear and channel-diagonal, so per-channel
    //! application is IDENTICAL to joint application -- splitting it is free.  An EXTRAPOLATOR is not: its
    //! coefficients come from a least-squares fit over a residual HISTORY, so independent fits per channel
    //! synthesise a state \f$(\sum c^\uparrow_i\rho^\uparrow_i,\ \sum c^\downarrow_i\rho^\downarrow_i)\f$ that
    //! never occurred on the trajectory -- each channel conserves its own charge, but the MOMENT becomes an
    //! arbitrary combination of history moments.  That is the MnO ejection (doc/SymmetryUpgradePlan.md §7
    //! step 7), and it is the same bug CLASS as the spin-blind ρ̃ mixer fixed at 041ddff3, one level up:
    //! that one was spin-blind in the FILTER, this one in the HISTORY.
    virtual tFieldExtrapolator* History() { return nullptr; }
};

//! One channel's contribution to a JOINT extrapolation: its own convergence gate, plus its block of the
//! shared error-overlap matrix.
struct StagedResidual
{
    double  resid=0.0;   //!< ‖out−in‖_∞ for THIS channel (the caller's convergence gate)
    rsmat_t B;           //!< \f$B^{(c)}_{ij}=\langle r^{(c)}_i,r^{(c)}_j\rangle\f$; EMPTY while priming / n<2
};

//! \brief The HISTORY face: a field mixer whose step is an EXTRAPOLATION over past iterations (Pulay/Broyden)
//! rather than a memoryless filter -- and which therefore splits its step into STAGING and APPLICATION so that
//! several channels can share ONE set of coefficients.
//!
//! The two-phase shape is forced by the physics: with two spin channels, ↑ cannot be mixed until ↓ has
//! contributed its block of B.  \c MixJointly below is the whole protocol; the single-channel \c MixField is
//! just its one-element case.
class tFieldExtrapolator : public virtual tFieldMixer
{
public:
    //! Phase 1: absorb \a out into this channel's history and return its block of the shared B (plus its
    //! residual).  The mix itself is DEFERRED -- nothing is committed until \c ApplyJoint.
    virtual StagedResidual StageResidual(const GField& out) = 0;
    //! Phase 2: commit the step with the coefficients \a c solved from the SUMMED B.  An EMPTY \a c means
    //! "no extrapolation this step" (priming, or history<2) -- the plain un-extrapolated step.
    virtual void ApplyJoint(const rvec_t& c) = 0;
    tFieldExtrapolator* History() override { return this; }
};

//! \brief ONE extrapolation over N channels: stage every channel, sum their blocks into a SINGLE bordered
//! solve, and apply the ONE coefficient vector to all of them.  Returns the WORST channel's residual.
//!
//! The joint inner product is simply the sum over channels, \f$B_{ij}=\sum_c\langle r^{(c)}_i,r^{(c)}_j
//! \rangle\f$ -- spin is just another irrep, exactly as the Fock-space DIIS already sums its B over spatial
//! irreps (SCFAcceleratorDIIS: "one B summed over both spins, one coefficient vector").  \c qchem.Math.DIIS
//! needs nothing new: it always documented that the CALLER owns the history and builds B with its own inner
//! product.
inline double MixJointly(const std::vector<tFieldExtrapolator*>& channels, const std::vector<GField>& fields)
{
    assert(channels.size()==fields.size() && !channels.empty());
    std::vector<StagedResidual> staged;
    double resid=0.0;
    for (size_t i=0;i<channels.size();++i)
    {
        staged.push_back(channels[i]->StageResidual(fields[i]));
        resid=std::max(resid,staged[i].resid);
    }
    const size_t n=staged[0].B.rows();
    rvec_t c;                                            // empty = no extrapolation (priming / n<2)
    if (n>=2)
    {
        rsmat_t B=blazem::zero<double>(n);
        for (const auto& s : staged)
        {
            assert(s.B.rows()==n && "MixJointly: the channels' histories are out of step");
            B+=s.B;
        }
        c=qchem::Math::DIIS::Coefficients(qchem::Math::DIIS::Bordered(B));
    }
    for (auto* ch : channels) ch->ApplyJoint(c);
    return resid;
}

// --- ΔG_Map arithmetic for Pulay (keys = integer G-index, consistent across iterations) ---
inline ΔG_Map MapSub(const ΔG_Map& a, const ΔG_Map& b)   // a - b
{
    ΔG_Map r=a;
    for (const auto& [k,v]:b) r[k]-=v;
    return r;
}
inline ΔG_Map MapCombine(const std::deque<ΔG_Map>& maps, const rvec_t& c)  // Σ cᵢ mapᵢ
{
    ΔG_Map r;
    for (size_t i=0;i<maps.size();++i)
        for (const auto& [k,v]:maps[i]) r[k]+=c[i]*v;
    return r;
}
inline ΔG_Map MapAdd(ΔG_Map a, const ΔG_Map& b)          // a + b (the ↑+↓ channel sum)
{
    for (const auto& [k,v]:b) a[k]+=v;
    return a;
}
inline ΔG_Map MapScale(ΔG_Map a, double s)               // s*a (the ½ of the (ρ±m)/2 channel rebuild)
{
    for (auto& [k,v]:a) v=s*v;
    return a;
}
//! s*(a + sign*b) for the raster shadows, empty unless BOTH arms answer on the same raster.
inline rvec_t RawCombine(const rvec_t& a, const rvec_t& b, double sign, double s)
{
    if (a.size()==0 || a.size()!=b.size()) return rvec_t{};
    rvec_t r=a; r+=sign*b; r*=s;
    return r;
}
inline double MapMaxAbs(const ΔG_Map& m)
{
    double x=0.0; for (const auto& [k,v]:m) x=std::max(x,std::abs(v)); return x;
}
inline double MapInnerRe(const ΔG_Map& a, const ΔG_Map& b)  // Re Σ_{G≠0} conj(aᵢ)·bⱼ (G=0 is never mixed)
{
    double s=0.0;
    for (const auto& [k,v]:a)
    {
        if (k==ivec3_t(0,0,0)) continue;
        auto it=b.find(k);
        if (it!=b.end()) s+=std::real(std::conj(v)*it->second);
    }
    return s;
}

} //namespace

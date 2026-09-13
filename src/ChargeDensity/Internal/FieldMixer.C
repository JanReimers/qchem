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
export import qchem.ChargeDensity.FourierDensity;         // FourierDensity, ΔG_Map + its field algebra (Projector3)
export import qchem.ChargeDensity.FourierMixCD;           // FourierMixCD (the field's PRESENTATION as a density)
export import qchem.ReciprocalLattice;                    // KerkerStep's |G|
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

//! The result of ONE Kerker step on a field: the mixed map, plus the two pieces of PROVENANCE the mixed
//! density's presentation carries for XC (see \c cDM_Sourced_CD): the N4 correction map (empty unless asked
//! for) and the realized mixing fraction \f$\alpha_{\rm eff}\f$.
struct KerkerStepResult
{
    ΔG_Map mix;              //!< \f$\tilde\rho_{in}+\alpha f_K(\tilde\rho_{out}-\tilde\rho_{in})\f$ over the union of both index sets
    ΔG_Map corr;             //!< N4: \f$\tilde\rho_{mix}-\tilde\rho_{out}\f$; empty when not formed
    double alphaEff=0.0;     //!< \f$\lVert\tilde\rho_{mix}-\tilde\rho_{in}\rVert_2/\lVert\tilde\rho_{out}-\tilde\rho_{in}\rVert_2\f$; 0 when converged
};

//! THE KERKER STEP, a pure function of fields: \f$\tilde\rho_{mix}(G)=\tilde\rho_{in}(G)+\alpha\,\frac{G^2}{G^2+G_0^2}
//! \,(\tilde\rho_{out}(G)-\tilde\rho_{in}(G))\f$.  \a G0 the Kerker screening wavevector (a.u.\f$^{-1}\f$;
//! \f$G_0\to0\f$ recovers plain linear mixing).  G=0 is mixed at full α (our ρ̃ is a fit-basis PROJECTION whose
//! (0,0,0) coefficient is shape-dependent, not the fixed \f$N/\Omega\f$; freezing it strands the XC's mean
//! density at the seed -- CP2K does the same).  \a withCuspCorrection (N4) also forms \c corr; FALSE leaves the
//! step bit-identical to its pre-N4 self.  Carries the GPW_KERKER_SPECTRUM residual-spectrum instrument.
KerkerStepResult KerkerStep(const ΔG_Map& in, const ΔG_Map& out, double alpha, double G0,
                            const ReciprocalLattice& recip, bool withCuspCorrection=false);

//! Build the PRESENTATION of a mixed field: the FourierMixCD the Fock is driven from, constructed whole from
//! the field the mixer owns and the step's provenance (never edited afterwards).  The N4 correction, when
//! formed, is presented as a field of ~ZERO net charge -- it is a difference of two densities of the same N.
inline std::shared_ptr<FourierMixCD> Present(const GField& f, const ReciprocalLattice& recip, double charge,
                                             const ΔG_Map& corr=ΔG_Map{}, double alphaEff=0.0)
{
    FourierMixCD::Extras x;
    x.raster=f.raster;
    if (!corr.empty()) x.xcCorrection=std::make_shared<const FourierMixCD>(corr, recip, 0.0);
    x.alphaEff=alphaEff;
    return std::make_shared<FourierMixCD>(f.tilde, recip, charge, std::move(x));
}

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
    //! The running mixed FIELD -- the mixer's own state (a caller recombining spin channels reads it here).
    virtual const GField& Field() const = 0;
    //! ...and its PRESENTATION as a density (what \c FockDensity returns).
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

//! s*(a + sign*b) for the raster shadows, empty unless BOTH arms answer on the same raster.
inline rvec_t RawCombine(const rvec_t& a, const rvec_t& b, double sign, double s)
{
    if (a.size()==0 || a.size()!=b.size()) return rvec_t{};
    rvec_t r=a; r+=sign*b; r*=s;
    return r;
}
//! The Pulay RESIDUAL metric \f$\mathrm{Re}\sum_{G\neq0}\overline{a_G}\,b_G\f$ -- G=0 is excluded because it is
//! never part of the mix's residual (the (0,0,0) coefficient is the charge, fixed by construction).  A
//! mixer-specific inner product, which is why it is here and not with the field algebra in Projector3.
inline double ResidualInnerRe(const ΔG_Map& a, const ΔG_Map& b)
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

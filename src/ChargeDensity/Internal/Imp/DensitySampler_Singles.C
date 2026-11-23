// File: ChargeDensity/Internal/Imp/DensitySampler_Singles.C  the SINGLES strategy: rho and H_xc both
// contracted through a cached table of SINGLE basis functions Phi_gi = chi_i(r_g).
//
// One implementation unit of module qchem.ChargeDensity.DensitySampler (extracted from
// qchem.Hamiltonian.Internal.PWTerms 2026-09-08 -- see that module's header for why).  Split 2026-09-08 out of a
// single 1213-line Imp/PWTerms.C (user: "PWTerms.C is huge, again doing too many things") into the
// interface-plus-many-Imp-units shape Internal/Terms.C has always had.  Helpers shared by more than
// one unit (NarrowExact, SampledField) live in the module INTERFACE's non-exported section, which is
// exactly what module-internal linkage is for -- they are visible to every unit of this module and to
// nothing outside it.
module;
#include <algorithm>   // std::min (the threaded quadrature's output-column blocking)
#include <cassert>
#include <cmath>       // std::abs (GroupResidualOf)
#include <complex>
#include <cstdlib>
#include <exception>   // std::exception_ptr (throw containment across the threaded Phi build)
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <type_traits>
#include <mutex>     // FitFor: the field-fit pool is shared by concurrent Fock blocks
#include <optional>    // the conditionally-charged sub-buckets of the H_xc quadrature
#include <stdexcept>
module qchem.ChargeDensity.Internal.DensitySampler;
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.RunPolicy;   // theRunPolicy().XCFromDM() -- the declared XC-feed deviation (N5)
import qchem.ChargeDensity;
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.G_FieldEvaluator;    // G_RasterTransform: the fit basis's FFT pair (RhoOnGrid, the BALL route)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.BasisSet.PointSetPullback;     // the block face the Becke mesh's pullback asks (READY 4)
import qchem.BasisSet.Transition_DFT_IBS;    // the (k+q, k) pair faces the transition Sample/Matrix go through (R3)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the cVxc_DFT engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Blaze;


namespace qchem::ChargeDensity
{

namespace {
// R2.1 (doc/SCFRobustnessPlan.md P2): max|f - Pf| over rho (and m), relative to max|rho| -- the PROJECTED rho \a scale, which a
// projector cannot inflate.  Empty m vectors = an unpolarized raster.
double GroupResidualOf(const rvec_t& rho0, const rvec_t& rho, const rvec_t& m0, const rvec_t& m, const rvec_t& scale)
{
    double d=0.0, r=0.0;
    for (size_t g=0; g<rho.size(); g++) d=std::max(d, std::abs(rho0[g]-rho[g]));
    for (size_t g=0; g<m.size();   g++) d=std::max(d, std::abs(m0[g]-m[g]));
    for (size_t g=0; g<scale.size(); g++) r=std::max(r, std::abs(scale[g]));
    return r>0.0 ? d/r : d;
}
} // namespace

// ---- SinglesDensitySampler: the pair-shared mesh + Phi tables + per-serial rho ----------------------------------

template <class T> tSinglesDensitySampler<T>::tSinglesDensitySampler(fit_t fit, BasisSet::FitQuadrature quad)
    : itsFit(std::move(fit))
    , itsQuad(std::move(quad))
{
    // A delta basis, so the delta fitter, by type (V1.34).  A named lvalue because the factory takes its argument by non-const
    // reference, as its siblings do.
    std::shared_ptr<const BasisSet::FIT_SF_ABS<T>> fb=itsFit;
    itsScalarFitter=Fitting::MakeDeltaScalarFitter(fb);
    assert(itsFit && "SinglesDensitySampler: the delta fit basis IS the quadrature -- it cannot be null");
    // The bundle's own invariants (fold partitions the mesh, sigmas need a fold, flags cover it) are
    // checked where the bundle becomes an object -- in the basis's ctor -- and the star-average announces
    // itself there too (providers self-report).  What IS this object's business is that the bundle it was
    // handed is the same one the basis holds: same length, so the fold's orbit indices and my coefficient
    // vectors index the same functions.
    assert((!itsQuad.GetMesh() || itsQuad.GetMesh()->size()==itsFit->GetNumFunctions()) &&
           "SinglesDensitySampler: the injected quadrature and the delta fit basis must be the SAME object "
           "-- one mesh point per fit function");
}

// THE FITTER'S PROJECTION FACE.  Since 2026-08-24 rho does NOT go to the fit basis: a projection is
// qcFitting's operation, and this strategy holds a fitter that already owns the fit basis's Phi handles
// for the adjoint direction (user).  So both directions of the SAME tensor now come off ONE object, which
// is the invariant this class's own header argues for one level up.
// CAPABILITY, not T (2026-11-12, lane H): the δ fitter's forward face, asked of the fitter itself.  An abstract -> abstract
// cross-cast on the real axis (the hub face \c DeltaScalarFit<double> does not name \c PairProjector<double>; the concrete
// fitter does), an up-cast on the complex one.  The real axis returned null here until then, so a molecule sampled ρ POINTWISE:
// every basis function at every mesh point on every sampling, 45 % of an LDA single point (OpenWork §4).
template <class T> const BasisSet::PairProjector<T>* tSinglesDensitySampler<T>::MatrixProjector() const
{
    return dynamic_cast<const BasisSet::PairProjector<T>*>(itsScalarFitter.get());
}
template <class T> const Fitting::FieldProjector& tSinglesDensitySampler<T>::FieldProjection() const
{
    return *itsScalarFitter;
}
template <class T> rvec_t tSinglesDensitySampler<T>::ProjectDM(const tDM_CD<T>& dm) const
{
    if (const BasisSet::PairProjector<T>* sp=MatrixProjector()) return dm.Project(*sp).Vector();   // the paired sampler: a delta fit's values
    return FieldProjection().Project(dm);
}

// The quadrature questions go THROUGH the basis, and since 2026-08-23 in its own FUNCTION vocabulary:
// f is an expansion over the delta basis (c_g = f(r_g), which is what makes a pointwise functional
// applicable to it at all), so its integral is the coefficients dotted with the functions' own integrals,
// Integral(sum_g c_g delta_g) = sum_g c_g <delta_g|1> = sum_g c_g w_g.  Algebraically the old
// Integrate(values) and, written as this loop, the SAME summation order as qcMesh::Integrate was -- which
// is what keeps the pinned energies bit-unmoved.
template <class T> const rvec_t& tSinglesDensitySampler<T>::FunctionIntegrals() const
{
    if (itsIntegrals.size()==0)
    {
        const vec_t<T> I=itsFit->Charge();
        itsIntegrals.resize(I.size());
        for (size_t a=0; a<I.size(); a++) itsIntegrals[a]=std::real(I[a]);
    }
    return itsIntegrals;
}
template <class T> double tSinglesDensitySampler<T>::Integrate(const rvec_t& f) const
{
    const rvec_t& I=FunctionIntegrals();
    assert(f.size()==I.size() && "SinglesDensitySampler::Integrate: one coefficient per fit function");
    double s=0.0;
    for (size_t a=0; a<f.size(); a++) s+=I[a]*f[a];
    return s;
}
template <class T> size_t tSinglesDensitySampler<T>::NumPoints() const {return itsFit->GetNumFunctions();}

// THE XC DM-rho REPAIR (doc/OpenWork.md, the factored-rho section).  Under rho-tilde mixing the density
// driving the Fock build is a G-space FIELD, so XC has to inverse-transform a truncated series at every
// mesh point -- 5.19 s per sampling on MnO (51% of the run) AND 8.5-18% of the points come back with
// rho<0, where the functionals guard `if (rho>0)` and silently contribute nothing.  But the mixer BUILT
// that field from a density matrix and now retains it (cDM_Sourced_CD), so the exact density is one
// cross-cast away: 0.042 s and rho>=0 by construction.  Hartree keeps the preconditioned field -- Poisson
// is LINEAR and diagonal in G -- while XC, a NONLINEAR POINTWISE functional, gets the cusp.  At the fixed
// point they agree, so this changes the SCF TRAJECTORY, not the answer.
// QCHEM_XC_DM_SOURCE=1 to arm it.  OPT-IN: a trajectory change must earn its place against banked recipes.
bool UseDMSource() {return theRunPolicy().XCFromDM();}
//! The exact density behind a field-backed one + THE DAMPING THAT FIELD APPLIED -- always together, because
//! taking the first without the second is the half-damped map (see cDM_Sourced_CD::EffectiveAlpha).
template <class T> struct ExactSource
{
    const tDM_CD<T>* cd=nullptr; double alpha=0.0;
    std::shared_ptr<const tChargeDensity<T>> corr;   //!< N4: rho_mix - rho[D] as a field (null = wholesale route)
    explicit operator bool() const {return cd;}
};
// THE CUSP-DEFICIT ROUTE (doc/OpenWork.md N4):
//     rho_XC = rho[D]_exact + IFT[rho_mix - rho[D]]
// instead of the WHOLESALE replacement (damped rho[D]) QCHEM_XC_DM_SOURCE takes.  The difference is not
// cosmetic: this form's BAND-LIMITED content is identically rho_mix -- the array Hartree uses -- so Kerker's
// f(G) reaches XC unmodified and there is NO alpha_eff to choose.  A flat alpha_eff cannot reproduce that
// SHAPE, and measurement (2026-08-25) says destroying the selectivity costs the MnO magnetic basin
// (-61.403 -> -45.529, Eee 13.5 -> 29.0 Ha).
// NOT A FLAG HERE, BY DESIGN (user, 2026-08-25): whether the correction exists is decided when the MIXER is
// BUILT (SCFParams::XCCuspDeficit -> KerkerMixerFactory), so CP2K parity is a property of the constructed
// object and the plain Kerker mixer stays bit-identical.  This code just uses what it was handed: a
// non-null XCCorrection() means the run asked for the correction route.
template <class T> ExactSource<T> ExactSourceOf(const qchem::ChargeDensity::tChargeDensity<T>* cd)
{
    if (!cd) return {};
    auto* src=dynamic_cast<const qchem::ChargeDensity::tDM_Sourced_CD<T>*>(cd);
    if (!src) return {};
    auto corr=src->XCCorrection();
    if (!corr && !UseDMSource()) return {};      // neither route armed
    return {src->DMSource().get(), src->EffectiveAlpha(), std::move(corr)};
}
// policy.xcDMMix overrides alpha_eff for CONTROLS only (=1 reproduces the undamped route); unset = use the
// mix's own.  policy.xcDMBoost scales it: alpha_eff came out ~0.20 on NaF and ~0.35 on MnO -- measured, but
// low against fractions those cells tolerate, and MnO converged 53 -> 39 iterations at boost 2.  NB f_K<=1
// forces alpha_eff<=alpha, so a boost >1/mean(f) puts XC ABOVE the mixer's own alpha -- defensible (XC's
// response kernel is finite at G->0, unlike Hartree's 4pi/G^2) but outside the preconditioner's bracket,
// hence a knob and not a default.  (Typed on RunPolicySpec since D-ENV 6c; the env variables are retired.)
double DMSourceMixOverride() { return qchem::theRunPolicy().XCDMMixOverride(); }
double DMSourceMixBoost()    { return qchem::theRunPolicy().XCDMBoost(); }
//! Blend \a fresh into \a running at the mix's own alpha (outside (0,1) => passthrough, which is the right
//! bootstrap AND the right answer for an unmixed or overshooting step).  Non-negativity survives: a convex
//! combination of non-negative rasters is non-negative.
void DampXCChannel(rvec_t& running, const rvec_t& fresh, double alphaEff)
{
    const double ov=DMSourceMixOverride();
    const double a =(ov>=0.0) ? ov : DMSourceMixBoost()*alphaEff;
    static const bool trace=qchem::Diagnostics::Enabled("xc_alpha");
    if (trace) std::cout<<"[XC alpha] alpha_eff="<<alphaEff<<(ov>=0.0?"  (OVERRIDDEN by policy.xcDMMix)":"")
                        <<"  boost="<<DMSourceMixBoost()
                        <<"  applied="<<((a>0.0&&a<1.0)?a:1.0)<<std::endl;
    if (a<=0.0 || a>=1.0 || running.size()!=fresh.size()) { running=fresh; return; }
    for (size_t g=0; g<fresh.size(); ++g) running[g]=(1.0-a)*running[g]+a*fresh[g];
}

// GPW_RHO_NEGATIVE=1: IS THE SAMPLED rho EVER NEGATIVE, AND WHAT DOES THAT COST?  A direct A/B between the
// two routes feeding this mesh, because they differ in KIND: the DM route gives rho = ||L^dag Phi||^2, a
// sum of squares, so rho>=0 BY CONSTRUCTION; the rho-tilde route inverse-transforms a TRUNCATED Fourier
// series, which RINGS, and rings NEGATIVE where the true density is sharpest -- the nuclear cusps this mesh
// exists to integrate.  And it is not loud: SlaterExchange::GetVxc guards `if (ro > 0.0)` and returns 0, so
// such a point contributes NOTHING to v_xc or E_xc with no diagnostic.  Hence the WEIGHTED MASS, not just
// the count -- that is the E_xc being silently dropped.
// Takes the QUADRATURE, not its mesh: the two masses are integrals, so it asks for them (R1.0) and never
// touches a weight -- which is what let the last Mesh() use out of this diagnostic.  It asks the STRATEGY
// rather than the fit basis (2026-08-23): integrating an expansion is the quadrature's question, and the
// basis now answers only the per-function pieces it is built from.
//! Does \a cd carry a retained-D source (QCHEM_XC_DM_SOURCE / the N4 cusp deficit)?  A bool-returning
//! wrapper so a caller declared ABOVE ExactSource's definition can still ask.
bool HasExactSource(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd) {return bool(ExactSourceOf(cd));}
//! ...and is it the N4 cusp-deficit route (a correction map), as opposed to the wholesale DM source?
bool HasXCCorrection(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd) {return bool(ExactSourceOf(cd).corr);}

//! The engine notes which ρ source it served XC from (the free helper below reports the sign census).
template <class T> void tSinglesDensitySampler<T>::NoteRoute(const char* route, XCRhoSource src) const
{
    itsSource=src;   // before the de-dup: a repeat route is still THIS sampling's source
    if (!itsRoutes.empty() && itsRoutes.back()==route) return;
    itsRoutes.push_back(route);
    // Contemporaneous, its own activity, inert without an open run report (the reporting rule).
    qchem::report::EmitAt("scf", "xcRhoRoutes", qchem::report::json(itsRoutes), qchem::report::Detail::Verbose);   // the json always records it (the Si gate reads it); the console only at Verbose -- it re-renders on every Fock build
}

void ReportNegativeRho(const DensitySamplerBase& q, const rvec_t& rho, const char* route)
{
    static const bool on=qchem::Diagnostics::Enabled("rho_negative");
    if (!on || rho.size()==0) return;
    rvec_t negOnly(rho.size(), 0.0), absRho(rho.size());
    size_t cnt=0; double minRho=0.0;
    for (size_t g=0; g<rho.size(); g++)
    {
        absRho[g]=std::fabs(rho[g]);
        if (rho[g]>=0.0) continue;
        cnt++; negOnly[g]=rho[g]; minRho=std::min(minRho,rho[g]);
    }
    const double negMass=q.Integrate(negOnly), absMass=q.Integrate(absRho);
    std::cout<<"[rho<0] route="<<route<<"  points="<<cnt<<" of "<<rho.size()
             <<" ("<<100.0*double(cnt)/double(rho.size())<<"%)"
             <<"  negative mass="<<negMass<<" e ("<<100.0*std::fabs(negMass)/std::max(absMass,1e-30)
             <<"% of integral|rho|)  min rho="<<minRho
             <<"   [these points contribute ZERO to E_xc -- GetVxc guards rho>0]"<<std::endl;
}

// rho at the mesh points, once per density serial for the WHOLE pair: the density GEMMs the cached
// tables against its private D (ProjectOnto; blocks not yet tabled self-evaluate pointwise -- first
// pass only).  A non-DM density (no DM face) falls back to the pointwise ScalarFunction sweep.
template <class T> const rvec_t& tSinglesDensitySampler<T>::Rho(const tChargeDensity<T>* cd) const
{
    assert(cd);
    // R2.9(i): the scalar and spin-resolved caches do not cross-invalidate (see the \warning on the
    // members).  An engine belongs to ONE xc/correlation pair, and a pair is either polarized or not, so
    // only one of the two routes is ever driven.  Pin it here rather than trusting the comment.
    assert(itsPolVersion==size_t(-1) && "SinglesDensitySampler: this engine already served RhoPol -- the scalar "
           "and spin-resolved rho caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsRhoVersion) return itsRho;
    itsRhoVersion=cd->Version();
    qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
    if (auto dm=dynamic_cast<const tDM_CD<T>*>(cd))
    {
        itsRho=ProjectDM(*dm);   // the density contracts its D into the fitter's handles
        NoteRoute("DM", XCRhoSource::DM);
        ReportNegativeRho(*this, itsRho, "DM");
    }
    else if (auto ex=ExactSourceOf(cd))
    {   // THE REPAIR, unpolarized sibling of the RhoPol branch below: the field retains the D it was mixed
        // from, so XC samples that instead of inverse-transforming a truncated series.
        if (ex.cd->Version()==itsSrcVersion)
            std::cerr<<"[XC DM-source] ** STALE: the mixed field advanced to serial "<<itsRhoVersion
                     <<" but its retained density matrix did NOT (still "<<itsSrcVersion
                     <<") -- V_xc is being built from the PREVIOUS iteration's D."<<std::endl;
        itsSrcVersion=ex.cd->Version();
        if (ex.corr)
        {   // N4: exact cusps + the band-limited difference.  No damping decision exists here by design.
            itsRho=ProjectDM(*ex.cd);
            const rvec_t c=FieldProjection().Project(*ex.corr);
            assert(c.size()==itsRho.size() && "XC cusp-deficit: correction and rho must share the mesh");
            for (size_t g=0; g<itsRho.size(); ++g) itsRho[g]+=c[g];
            NoteRoute("cusp-deficit", XCRhoSource::CuspDeficit);
            ReportNegativeRho(*this, itsRho, "cusp-deficit");
        }
        else
        {
            DampXCChannel(itsXCMix, ProjectDM(*ex.cd), ex.alpha);
            itsRho=itsXCMix;   // the running mix lives in its OWN buffer -- see the \warning on itsXCMix
            NoteRoute("DM-source", XCRhoSource::DMSource);
            ReportNegativeRho(*this, itsRho, "DM-source");
        }
    }
    else
    {
        itsRho=FieldProjection().Project(*cd);   // non-DM (mixed rho-tilde / seed): the fitter projects the FIELD
        NoteRoute("matrix-free", XCRhoSource::Fitted);
        ReportNegativeRho(*this, itsRho, "matrix-free");
    }
    {
        const rvec_t raw=itsRho;      // R2.1: what XC was handed, measured against the projector it is about to apply
        itsQuad.Symmetrize(itsRho);   // §6a W1: the orbit-mean projector, straight off the FoldedMesh I hold
                                      // (no-op on a free run, so no caller asks whether symmetry was imposed)
        itsGroupResidual=GroupResidualOf(raw, itsRho, {}, {}, itsRho);
    }
    return itsRho;
}

// δρ of a transition density at my points (R2; pair-shaped since R3 step 3): each block pair through the δ basis's
// Φ tables (Transition_Overlap3C) -- the same tables as Rho, deliberately WITHOUT itsQuad.Symmetrize (see the face)
// and without touching the caches.
namespace {
//! The δ basis's transition pair as the density's projector (q unused: the point route samples δρ itself).
class DeltaTransitionProjector : public TransitionProjector
{
public:
    DeltaTransitionProjector(const BasisSet::Transition_Overlap3C& f, size_t n) : itsF(f), itsN(n) {}
    virtual cvec_t Forward(const cobs_t& bra, const cobs_t& ket, const mat_t<dcmplx>& dD, const rvec3_t&) const override
    {
        return itsF.TransitionForward(dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(bra),
                                      dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(ket), dD);
    }
    virtual size_t NumCoefficients() const override {return itsN;}
private:
    const BasisSet::Transition_Overlap3C& itsF;
    size_t                                itsN;
};
const BasisSet::Transition_Overlap3C& TransitionFace(const BasisSet::cFIT_SF_ABS& fit)
{
    auto* f=dynamic_cast<const BasisSet::Transition_Overlap3C*>(&fit);   // abstract -> abstract
    if (!f) throw std::logic_error("SinglesDensitySampler: this point fit basis has no Transition_Overlap3C face");
    return *f;
}
} // namespace

cvec_t SinglesDensitySampler::Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const
{
    const auto* ch=delta.Channel(s);
    if (!ch) throw std::logic_error("SinglesDensitySampler::Sample: the transition density does not resolve this spin channel");
    auto* p=dynamic_cast<const ProjectableTransition*>(ch);
    if (!p) throw std::logic_error("SinglesDensitySampler::Sample: this transition density cannot be projected pair by pair");
    qchem::report::Timed timed("response: XC-mesh delta-rho sampling");
    return p->ProjectOnto(DeltaTransitionProjector(TransitionFace(*itsFit), NumPoints()));
}

mat_t<dcmplx> SinglesDensitySampler::Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t&) const
{
    return TransitionFace(*itsFit).TransitionAdjoint(dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(*bra),
                                                     dynamic_cast<const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&>(*ket), v);
}

// The real-block ensure siblings (3c-3): build the real block's OWN typed table first (PhiR), then the
// shared sampling path -- the exact mirror of the complex ensureBlock argument.

// The spin-resolved sibling of Rho: the {up,down} PAIR is cached under ONE density serial (a polarized
// density's Version() forwards to its Up child -- a single scalar cache would hand the Up raster to the
// Down channel).  A polarized composite answers per channel VIEW (each channel GEMMs its own D against the
// SHARED Phi tables); a spin-agnostic density (the seed) collapses to rho/2 per channel, so the first
// iterations run the exact unpolarized collapse (v^sigma(rho/2,rho/2)=v^P(rho)).
template <class T> const rvec_t& tSinglesDensitySampler<T>::RhoPol(const tChargeDensity<T>* cd, const Spin& s) const
{
    assert(cd);
    assert(itsRhoVersion==size_t(-1) && "SinglesDensitySampler: this engine already served the scalar Rho -- the "
           "two rho caches have no cross-invalidation, so one of them would go stale");
    assert(s!=Spin::None && "SinglesDensitySampler::RhoPol: ask for a channel, not the total");
    if (cd->Version()!=itsPolVersion)
    {
        itsPolVersion=cd->Version();
        qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
        // THE CHANNELS THROUGH THE FACE (V1.37): a polarized composite answers its Up/Down VIEWS, each a
        // tDM_CD<T> that GEMMs its own D against the SHARED Phi tables; the matrix-free pair below answers its
        // channel objects; a spin-agnostic density answers null and collapses to rho/2 per channel.
        const tChargeDensity<T>* up=ChannelOf(cd, Spin::Up  );
        const tChargeDensity<T>* dn=ChannelOf(cd, Spin::Down);
        if (auto* dmUp=dynamic_cast<const tDM_CD<T>*>(up))
        {
            auto* dmDn=dynamic_cast<const tDM_CD<T>*>(dn);
            assert(dmDn && "a D-backed Up channel beside a matrix-free Down channel");
            itsRhoUp=ProjectDM(*dmUp);
            itsRhoDn=ProjectDM(*dmDn);
            NoteRoute("DM(up)", XCRhoSource::DM);
            ReportNegativeRho(*this, itsRhoUp, "DM(up)");
            NoteRoute("DM(dn)", XCRhoSource::DM);
            ReportNegativeRho(*this, itsRhoDn, "DM(dn)");
        }
        else if (up && dn)
        {   // MATRIX-FREE spin-resolved density: the seed (PolarizedSeedCD, SCFSeedingPlan §10) at
            // iteration 0, and -- the expensive case -- the ρ̃-MIXED density (PolarizedMixCD over
            // FourierMixCD) on EVERY Kerker/Pulay iteration.  Neither carries a D, so both batch through
            // the plain tChargeDensity face instead of the Phi GEMM: an image-summed atomic sum for the
            // seed, a batched inverse FT over the whole {G} for the mixer.  Its OWN bucket because it is
            // a different algorithm at a different cadence -- lumping it into the GEMM hid the fact that
            // the mixed-density sampling, not the GEMM, was the iteration's largest XC cost.
            // THE REPAIR: each channel of a rho-tilde-mixed density retains the DM-backed channel it was
            // mixed from, so ask for the exact one FIRST and fall back per CHANNEL -- the seed has no source
            // on either while a mixed density has one on both, so a pair-level test would be right today and
            // wrong the first time they differ.
            const ExactSource<T> exUp=ExactSourceOf(up);
            const ExactSource<T> exDn=ExactSourceOf(dn);
            if (exUp && exDn)
            {
                // STALENESS GUARD (see itsSrcVersion): the field's serial just advanced, so the SOURCE's
                // must have too, or XC gets last iteration's D under a cache that believes it is fresh.
                if (exUp.cd->Version()==itsSrcVersion)
                    std::cerr<<"[XC DM-source] ** STALE: the mixed field advanced to serial "<<itsPolVersion
                             <<" but its retained density matrix did NOT (still "<<itsSrcVersion
                             <<") -- V_xc is being built from the PREVIOUS iteration's D."<<std::endl;
                itsSrcVersion=exUp.cd->Version();
                rvec_t up=ProjectDM(*exUp.cd);
                rvec_t dn=ProjectDM(*exDn.cd);
                if (exUp.corr && exDn.corr)
                {   // N4, per channel: each channel's own rho_mix - rho[D] (the polarized split gives the
                    // mixer one FourierMixCD per spin, so the corrections are already channel-private).
                    const rvec_t cu=FieldProjection().Project(*exUp.corr), cd_=FieldProjection().Project(*exDn.corr);
                    assert(cu.size()==up.size() && cd_.size()==dn.size());
                    for (size_t g=0; g<up.size(); ++g) {up[g]+=cu[g]; dn[g]+=cd_[g];}
                    itsRhoUp=up; itsRhoDn=dn;
                    NoteRoute("cusp-deficit(up)", XCRhoSource::CuspDeficit);
                    ReportNegativeRho(*this, itsRhoUp, "cusp-deficit(up)");
                    NoteRoute("cusp-deficit(dn)", XCRhoSource::CuspDeficit);
                    ReportNegativeRho(*this, itsRhoDn, "cusp-deficit(dn)");
                }
                else
                {
                    DampXCChannel(itsXCMixUp, up, exUp.alpha);   // match the damping Hartree gets, so the map
                    DampXCChannel(itsXCMixDn, dn, exDn.alpha);   //   is not half-damped
                    itsRhoUp=itsXCMixUp; itsRhoDn=itsXCMixDn;
                    NoteRoute("DM-source(up)", XCRhoSource::DMSource);
                    ReportNegativeRho(*this, itsRhoUp, "DM-source(up)");
                    NoteRoute("DM-source(dn)", XCRhoSource::DMSource);
                    ReportNegativeRho(*this, itsRhoDn, "DM-source(dn)");
                }
            }
            else
            {
            qchem::report::Timed seed("scf: XC-mesh rho sampling (matrix-free density)");
            itsRhoUp=FieldProjection().Project(*up);
            itsRhoDn=FieldProjection().Project(*dn);
            NoteRoute("matrix-free(up)", XCRhoSource::Fitted);
            ReportNegativeRho(*this, itsRhoUp, "matrix-free(up)");
            NoteRoute("matrix-free(dn)", XCRhoSource::Fitted);
            ReportNegativeRho(*this, itsRhoDn, "matrix-free(dn)");
            }
        }
        else
        {   // spin-agnostic seed: rho_up=rho_down=rho/2 (the molecular HalfDensity rule, cd85d13c)
            if (auto dm=dynamic_cast<const tDM_CD<T>*>(cd))
            {
                itsRhoUp=ProjectDM(*dm);
                itsSource=XCRhoSource::DM;
            }
            else
            {
                itsRhoUp=FieldProjection().Project(*cd);
                itsSource=XCRhoSource::Fitted;
            }
            itsRhoUp*=0.5;
            itsRhoDn=itsRhoUp;
        }
        {
            // Project the (ρ,m) PAIR over the injected quadrature: magnetic (σ tags -> ρ even, m odd with
            // the flip-fixed zeros), grey (each argument averaged independently) or free (a no-op).  The
            // three-way decision belongs to the FoldedMesh, which owns the fold, the σ tags and the
            // flip-fixed audit -- this class only knows WHICH pair to hand it.
            rvec_t rho = itsRhoUp + itsRhoDn;
            rvec_t m   = itsRhoUp - itsRhoDn;
            const rvec_t rho0=rho, m0=m;   // R2.1: the pair XC was handed, measured against the projector it is about to apply
            itsQuad.SymmetrizeSpin(rho, m);
            itsGroupResidual=GroupResidualOf(rho0, rho, m0, m, rho);
            itsRhoUp = 0.5*(rho + m);
            itsRhoDn = 0.5*(rho - m);
        }
    }
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

// The ONE place the injected partition is read: Integral w_A f over each site block.  Not a fit-basis
// question and no longer asked of one -- the mesh arrives from the factory that built it for the basis,
// so this and the basis index the SAME points in the SAME order (one object, two collaborators).
template <class T> rvec_t tSinglesDensitySampler<T>::SiteIntegrals(const rvec_t& f) const
{
    if (!itsQuad.GetMesh() || itsQuad.GetMesh()->NSites()==0)
    {   // No site partition on this mesh -- legitimate for a uniform grid, a DEFECT for an atom-centred one.  The facade's
        // result line says which (SolidCalculation::Converge), so the sampler stays silent: it only answers "none".
        return rvec_t();
    }
    assert(f.size()==itsQuad.GetMesh()->size() && "SinglesDensitySampler: the injected quadrature and the fit "
           "basis must be the same object -- one field value per mesh point");
    return qcMesh::SiteIntegrals(*itsQuad.GetMesh(), f);
}

// <i|v|j>: ASK THE BASIS.  It owns the points, the weights and the Phi table, so the whole quadrature
// -- Phi^dag diag(w v) Phi -- is its operation; this strategy only decides WHICH v to hand it.  That is
// what closed the last weight/coordinate escape (doc/CleanupCandidates.md R1.0 increment 2).
// <i|v|j> THROUGH THE FITTER, exactly as the molecular XC term does it (the Liskov conformance,
// doc/CleanupCandidates.md R1.0): fit the sampled field, then contract against this block.  The fit is
// re-done only when v CHANGES -- the Fock build calls this once per block with the same v, and DoFit
// on a delta basis is a copy, but re-copying per block would still be per-block work for nothing.  One fit PER v
// (FitFor), so the up and down blocks of a concurrent Fock loop never refit each other's fitter mid-contraction.
template <class T> std::shared_ptr<const Fitting::DeltaScalarFit<T>> tSinglesDensitySampler<T>::FitFor(const rvec_t& v) const
{
    auto same=[&](const rvec_t& a)
    {
        if (a.size()!=v.size()) return false;
        for (size_t g=0; g<v.size(); g++) if (a[g]!=v[g]) return false;
        return true;
    };
    std::lock_guard<std::mutex> lock(itsFieldFitsLock.m);
    for (const auto& f : itsFieldFits) if (same(f.v)) return f.fitter;
    std::shared_ptr<Fitting::DeltaScalarFit<T>> fitter;
    if (itsFieldFits.size()>=kFieldFits)                       // recycle the OLDEST slot no block is reading
        for (auto i=itsFieldFits.begin(); i!=itsFieldFits.end(); ++i)
            if (i->fitter.use_count()==1) {fitter=std::move(i->fitter); itsFieldFits.erase(i); break;}
    if (!fitter)
    {
        std::shared_ptr<const BasisSet::FIT_SF_ABS<T>> fb=itsFit;   // Factory takes a non-const lvalue (see the ctor)
        fitter=Fitting::MakeDeltaScalarFitter(fb);
    }
    fitter->DoFit(SampledField(v, NumPoints()));                // under the lock: on a delta basis the fit is a copy
    itsFieldFits.push_back({v, fitter});
    return fitter;
}

template <class T> template <class U> hmat_t<U> tSinglesDensitySampler<T>::MatrixT(const tobs_t<U>* bs, const rvec_t& v) const
{
    const std::shared_ptr<const Fitting::DeltaScalarFit<T>> fitter=FitFor(v);   // held: its slot cannot be recycled
    // TWO scalars on the face, and they DIFFER here (3c-3): the block is U, but the delta fit basis this
    // strategy holds is the run's periodic one, so the fit axis is dcmplx for both.  <U> alone would name
    // <U,U>, which for a real TRIM block is a face no fitter in the tree declares.  And the DFT cross-cast
    // is the XC strategy's, for the reason given on FitContraction: this layer is DFT by construction.
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,T>&>(*bs);
    const Fitting::FitContraction<U,T>& contraction=*fitter;   // a block scalar this fit lacks does not compile (V1.34)
    return contraction.Overlap(orb);
}
template <class T> hmat_t<T> tSinglesDensitySampler<T>::Matrix(const tobs_t<T>* bs, const rvec_t& v) const {return MatrixT<T>(bs,v);}
rsmat_t SinglesDensitySampler::Matrix(const robs_t* bs, const rvec_t& v) const {return MatrixT<double>(bs,v);}

// ---- THE DERIVATIVES ON A FREE ATOM-CENTRED MESH (lane B 2026-11-19, READY 4; the controller-approved split of the design note) ------------
// E_xc = Σ_g w_g e(ρ_g), w_g = w_q W(R_owner + u_g).  At fixed D, dE = Σ_g ∂w_g e_g + Σ_g w_g v_g ∂ρ_g: the second sum per block (the
// functions moving, the points riding their owners, the weights HELD: MakeMatrixPullback / MakeMatrixVirial), the first the quadrature's own
// (MakeQuadraturePullback, the mesh builder's weight response; the XC term contracts it with e).
template <class T> std::string tSinglesDensitySampler<T>::MeshDerivativeGap() const
{
    const auto& m=itsQuad.GetMesh();
    if (!m) return "the XC quadrature has no mesh injected";
    if (itsQuad.HasFold()) return "an IMPOSED-symmetry XC mesh (orbit-rotated grids, a star-averaged rho) gives no geometry derivative yet";
    if (m->OwnerOffsets().size()!=m->size() || m->NSites()==0) return "the XC mesh records no point owners (not an atom-centred build)";
    if (!m->WeightResponse()) return "the XC mesh's builder installed no weight response";
    return "";
}
template <class T> bool tSinglesDensitySampler<T>::HasQuadratureResponse() const {return MeshDerivativeGap().empty();}
template <class T> Outcome<typename tDensitySampler<T>::QuadraturePullback,AD::NoGradient>
tSinglesDensitySampler<T>::MakeQuadraturePullback(const rvec_t& f, const Structure* st) const
{
    using R=Outcome<typename tDensitySampler<T>::QuadraturePullback,AD::NoGradient>;
    if (const std::string gap=MeshDerivativeGap(); !gap.empty()) return R::Fail({gap});
    const qcMesh::Mesh& m=*itsQuad.GetMesh();
    if (!st || st->GetNumAtoms()!=m.NSites()) return R::Fail({"the XC mesh's site blocks are not the structure's atoms"});
    const auto r=m.WeightResponse()->Pullback(m, f);
    return R::Ok({r.gradient, r.virial});
}

// THE BLOCK's FACE.  A cross-cast per call, abstract -> abstract, as Matrix's own DFT cross-cast: this sampler is built ONCE per run (by
// MakeDensitySampler, from the fit basis and its quadrature) but serves every block the term hands it per call, Bloch and real TRIM alike,
// so a capability handed at construction would have to be a per-block table keyed by the block -- the probing it would replace (the
// controller's condition (a), 2026-11-19).  A block without the face is a named refusal (CanT false, PointSetOf's Fail), never a null deref.
template <class T> template <class U> bool tSinglesDensitySampler<T>::CanT(const tobs_t<U>* bs) const
{
    return MeshDerivativeGap().empty() && dynamic_cast<const BasisSet::PointSetPullback<U>*>(bs)!=nullptr;
}
template <class T> template <class U> Outcome<typename BasisSet::PointSetPullback<U>::Result,AD::NoGradient>
tSinglesDensitySampler<T>::PointSetOf(const tobs_t<U>* bs, const rvec_t& v, const hmat_t<U>& D, const Structure* st) const
{
    using R=Outcome<typename BasisSet::PointSetPullback<U>::Result,AD::NoGradient>;
    if (const std::string gap=MeshDerivativeGap(); !gap.empty()) return R::Fail({gap});
    const auto* face=dynamic_cast<const BasisSet::PointSetPullback<U>*>(bs);
    if (!face) return R::Fail({"the orbital block "+bs->Name()+" gives no point-set pullback"});
    const qcMesh::Mesh& m=*itsQuad.GetMesh();
    if (v.size()!=m.size()) return R::Fail({"the XC field is not on the mesh"});
    if (!st || st->GetNumAtoms()!=m.NSites()) return R::Fail({"the XC mesh's site blocks are not the structure's atoms"});
    rvec_t c(m.size());
    const rvec_t& w=m.Weights();
    for (size_t g=0; g<c.size(); g++) c[g]=w[g]*v[g];
    return face->MakePointSetPullback(m.Points(), c, D, st);
}
template <class T> template <class P> rvec3vec_t tSinglesDensitySampler<T>::FinishPullback(const P& p) const
{
    rvec3vec_t g=p.gradient;
    const qcMesh::Mesh& m=*itsQuad.GetMesh();
    for (size_t a=0; a<m.NSites(); a++)
        for (size_t q=m.SiteBegin(a); q<m.SiteEnd(a); q++) g[a]+=p.perPoint[q];   // the point rides its owner
    return g;
}
template <class T> template <class P> rmat3d_t tSinglesDensitySampler<T>::FinishVirial(const P& p) const
{
    rmat3d_t V=p.moment;
    const rvec3vec_t& Uq=itsQuad.GetMesh()->OwnerOffsets();
    for (size_t q=0; q<Uq.size(); q++)   // the point is X_owner + u (u not strained): moment - perPoint (x) u
    {
        const double pp[3]={p.perPoint[q].x, p.perPoint[q].y, p.perPoint[q].z}, uu[3]={Uq[q].x, Uq[q].y, Uq[q].z};
        for (size_t a=0;a<3;a++) for (size_t b=0;b<3;b++) V(a+1,b+1)-=pp[a]*uu[b];
    }
    return V;
}
template <class T> template <class U> Outcome<rvec3vec_t,AD::NoGradient>
tSinglesDensitySampler<T>::PullbackT(const tobs_t<U>* bs, const rvec_t& v, const hmat_t<U>& D, const Structure* st) const
{
    const auto p=PointSetOf(bs, v, D, st);
    if (!p) return Outcome<rvec3vec_t,AD::NoGradient>::Fail(p.Error());
    return Outcome<rvec3vec_t,AD::NoGradient>::Ok(FinishPullback(*p));
}
template <class T> template <class U> Outcome<rmat3d_t,AD::NoGradient>
tSinglesDensitySampler<T>::VirialT(const tobs_t<U>* bs, const rvec_t& v, const hmat_t<U>& D, const Structure* st) const
{
    const auto p=PointSetOf(bs, v, D, st);
    if (!p) return Outcome<rmat3d_t,AD::NoGradient>::Fail(p.Error());
    return Outcome<rmat3d_t,AD::NoGradient>::Ok(FinishVirial(*p));
}

// THE k-FOLD of the whole-density calls (user ruling 2026-11-22 "yes").  Each block (either scalar) is one entry; entries whose blocks are
// KFoldable with one non-zero domain AND see the same field (one spin's v) form a group, and a group of two or more goes through its first
// block's MakePointSetPullbackFolded with every member's D widened onto the complex face.  The rest go alone, exactly as the per-block calls.
// Linear end to end: the finish (owners, the u moment) of a summed result is the sum of the finishes.
template <class T> template <class Res, class Finish, class Add> Outcome<Res,AD::NoGradient>
tSinglesDensitySampler<T>::FoldedBlocksT(const SampledBlocks<T>& b, const Structure* st, Res total, Finish finish, Add add) const
{
    using R=Outcome<Res,AD::NoGradient>;
    if (const std::string gap=MeshDerivativeGap(); !gap.empty()) return R::Fail({gap});
    struct Entry
    {
        const tobs_t<T>*      bs=nullptr;  const hmat_t<T>*      D=nullptr;   // exactly one of the two pairs is set
        const tobs_t<double>* rbs=nullptr; const hmat_t<double>* rD=nullptr;
        const rvec_t* v=nullptr;
        const BasisSet::KFoldable* kf=nullptr;
    };
    std::vector<Entry> e;
    for (const auto& x : b.blocks) e.push_back({x.bs, x.D, nullptr, nullptr, x.v, dynamic_cast<const BasisSet::KFoldable*>(x.bs)});
    for (const auto& x : b.real)   e.push_back({nullptr, nullptr, x.bs, x.D, x.v, dynamic_cast<const BasisSet::KFoldable*>(x.bs)});
    auto alone=[&](const Entry& x)->R
    {
        if (x.bs) { const auto p=PointSetOf(x.bs, *x.v, *x.D, st); if (!p) return R::Fail(p.Error()); return R::Ok(finish(*p)); }
        if constexpr (std::is_same_v<T,double>)
        { const auto p=PointSetOf(x.rbs, *x.v, *x.rD, st); if (!p) return R::Fail(p.Error()); return R::Ok(finish(*p)); }
        else
        { const auto p=PointSetOf<double>(x.rbs, *x.v, *x.rD, st); if (!p) return R::Fail(p.Error()); return R::Ok(finish(*p)); }
    };
    std::vector<bool> done(e.size(), false);
    for (size_t i=0;i<e.size();i++)
    {
        if (done[i]) continue;
        done[i]=true;
        std::vector<size_t> grp{i};
        if (e[i].kf && e[i].kf->FoldDomain()[0]!=0)
            for (size_t j=i+1;j<e.size();j++)
                if (!done[j] && e[j].kf && e[j].v==e[i].v && e[j].kf->FoldDomain()==e[i].kf->FoldDomain()) { grp.push_back(j); done[j]=true; }
        if (grp.size()==1) { const R r=alone(e[i]); if (!r) return r; add(total, *r); continue; }
        // The group: the lead's face, the weights c = w v, every member's D on the complex face.
        const qcMesh::Mesh& m=*itsQuad.GetMesh();
        const rvec_t& v=*e[i].v;
        if (v.size()!=m.size()) return R::Fail({"the XC field is not on the mesh"});
        if (!st || st->GetNumAtoms()!=m.NSites()) return R::Fail({"the XC mesh's site blocks are not the structure's atoms"});
        rvec_t c(m.size());
        const rvec_t& w=m.Weights();
        for (size_t g=0; g<c.size(); g++) c[g]=w[g]*v[g];
        std::vector<chmat_t> widened; widened.reserve(grp.size());   // reserved: the terms point into it
        std::vector<BasisSet::KFoldTerm> terms;
        for (size_t j : grp)
        {
            const Entry& x=e[j];
            if (x.bs)
            {
                if constexpr (std::is_same_v<T,dcmplx>) { terms.push_back({x.kf, x.D}); continue; }
                else { chmat_t d(x.D->rows()); for (size_t r=0;r<d.rows();r++) for (size_t s=r;s<d.rows();s++) d(r,s)=dcmplx((*x.D)(r,s),0.0);
                       widened.push_back(std::move(d)); }
            }
            else { chmat_t d(x.rD->rows()); for (size_t r=0;r<d.rows();r++) for (size_t s=r;s<d.rows();s++) d(r,s)=dcmplx((*x.rD)(r,s),0.0);
                   widened.push_back(std::move(d)); }
            terms.push_back({x.kf, &widened.back()});
        }
        auto viaLead=[&](const auto* face)->R
        {
            if (!face) return R::Fail({"the orbital block gives no point-set pullback"});
            const auto p=face->MakePointSetPullbackFolded(m.Points(), c, terms, st);
            if (!p) return R::Fail(p.Error());
            return R::Ok(finish(*p));
        };
        const R r= e[i].bs ? viaLead(dynamic_cast<const BasisSet::PointSetPullback<T>*>(e[i].bs))
                           : viaLead(dynamic_cast<const BasisSet::PointSetPullback<double>*>(e[i].rbs));
        if (!r) return r;
        add(total, *r);
    }
    return R::Ok(std::move(total));
}
template <class T> Outcome<rvec3vec_t,AD::NoGradient> tSinglesDensitySampler<T>::MakeMatrixPullbackBlocks(const SampledBlocks<T>& b, const Structure* st) const
{
    return FoldedBlocksT<rvec3vec_t>(b, st, rvec3vec_t{},
        [this](const auto& p){return FinishPullback(p);},
        [](rvec3vec_t& F, const rvec3vec_t& f)
        {
            if (F.size()==0) { F=f; return; }
            if (f.size()!=F.size()) throw std::logic_error("MakeMatrixPullbackBlocks: two blocks answered for different atom counts");
            for (size_t A=0;A<F.size();A++) F[A]+=f[A];
        });
}
template <class T> Outcome<rmat3d_t,AD::NoGradient> tSinglesDensitySampler<T>::MakeMatrixVirialBlocks(const SampledBlocks<T>& b, const Structure* st) const
{
    return FoldedBlocksT<rmat3d_t>(b, st, rmat3d_t(0,0,0, 0,0,0, 0,0,0),
        [this](const auto& p){return FinishVirial(p);},
        [](rmat3d_t& V, const rmat3d_t& v){V+=v;});
}
template <class T> bool tSinglesDensitySampler<T>::CanPullbackMatrix(const tobs_t<T>* bs) const {return CanT<T>(bs);}
template <class T> bool tSinglesDensitySampler<T>::CanVirialMatrix  (const tobs_t<T>* bs) const {return CanT<T>(bs);}
template <class T> Outcome<rvec3vec_t,AD::NoGradient> tSinglesDensitySampler<T>::MakeMatrixPullback(const tobs_t<T>* bs, const rvec_t& v, const hmat_t<T>& D, const Structure* st) const
{return PullbackT<T>(bs, v, D, st);}
template <class T> Outcome<rmat3d_t,AD::NoGradient> tSinglesDensitySampler<T>::MakeMatrixVirial(const tobs_t<T>* bs, const rvec_t& v, const hmat_t<T>& D, const Structure* st) const
{return VirialT<T>(bs, v, D, st);}
bool SinglesDensitySampler::CanPullbackMatrix(const robs_t* bs) const {return CanT<double>(bs);}
bool SinglesDensitySampler::CanVirialMatrix  (const robs_t* bs) const {return CanT<double>(bs);}
Outcome<rvec3vec_t,AD::NoGradient> SinglesDensitySampler::MakeMatrixPullback(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const
{return PullbackT<double>(bs, v, D, st);}
Outcome<rmat3d_t,AD::NoGradient> SinglesDensitySampler::MakeMatrixVirial(const robs_t* bs, const rvec_t& v, const hmat_t<double>& D, const Structure* st) const
{return VirialT<double>(bs, v, D, st);}

template class tSinglesDensitySampler<double>;   // the molecular / atom path (V-XC-ONE-TERM)
template class tSinglesDensitySampler<dcmplx>;   // the periodic path

// ---- THE DENSITY-GRADIENT STRATEGY (GGA stages 2-3, the point route) ------------------------------------------------------------------

template <class T> DeltaGradientProjector<T>::DeltaGradientProjector(std::shared_ptr<const BasisSet::FIT_SF_ABS<T>> fb, const qcMesh::Mesh& mesh)
    : itsFit(std::move(fb))
    , itsPoints(mesh.Points())
{
    if (!dynamic_cast<const BasisSet::Integrals_GradOverlap3C<T,T>*>(itsFit.get()))   // abstract -> abstract
        throw std::invalid_argument("DeltaGradientProjector: the fit basis "+itsFit->Name()+" has no gradient 3-centre overlap "
                                    "(BasisSet::Integrals_GradOverlap3C), so this route cannot sample the density gradient");
    if (itsPoints.size()!=itsFit->GetNumFunctions())
        throw std::logic_error("DeltaGradientProjector: the mesh and the delta basis must be one object -- one point per function");
    const rvec_t& w=mesh.Weights();
    const size_t n=w.size();
    rvec_t w3(3*n);
    for (size_t a=0; a<3; a++) for (size_t g=0; g<n; g++) w3[a*n+g]=w[g];
    itsW3=std::make_shared<const rvec_t>(std::move(w3));
}

template <class T> template <class U>
const ScreenedMatrixIntegrator<U>& DeltaGradientProjector<T>::Integrator(const BasisSet::Orbital_DFT_IBS<U,T>& orb) const
{
    std::lock_guard<std::mutex> lock(itsIntLock.m);   // filled from the parallel Fock pass
    return Cache<U>().GetOrFill(orb.Instance(), [&]
    {
        const auto* face=dynamic_cast<const BasisSet::Integrals_GradOverlap3C<U,T>*>(itsFit.get());   // abstract -> abstract
        if (!face) throw std::logic_error("DeltaGradientProjector: the fit basis has no gradient 3-centre overlap for this block scalar");
        return ScreenedMatrixIntegrator<U>(face->GradOverlap3C(orb), itsW3);
    });
}

template <class T> Projections DeltaGradientProjector<T>::ProjectField(const ScalarFunction<double>& f) const
{
    const size_t n=itsPoints.size();
    rvec_t c(3*n);
    for (size_t g=0; g<n; g++)
    {
        const rvec3_t d=f.Gradient(itsPoints[g]);
        c[g]=d.x; c[n+g]=d.y; c[2*n+g]=d.z;
    }
    return Projections(std::move(c));
}

template <class T, class Base> tSinglesGradientSampler<T,Base>::tSinglesGradientSampler(typename tSinglesDensitySampler<T>::fit_t fit,
                                                                                         BasisSet::FitQuadrature quad)
    : Base(fit, quad)
    , itsGrad(fit, *quad.GetMesh())
{
    if (this->itsQuad.HasFold())
        throw std::invalid_argument("tSinglesGradientSampler: a FOLDED quadrature star-averages rho over orbits, and the density gradient "
                                    "is a vector field the fold would have to rotate -- no fold does that yet, so a GGA runs free");
}

template <class T, class Base> rvec3vec_t tSinglesGradientSampler<T,Base>::Unpack(const rvec_t& c) const
{
    const size_t n=this->NumPoints();
    if (c.size()!=3*n) throw std::logic_error("tSinglesGradientSampler: a gradient projection is not three components per point");
    rvec3vec_t out(n);
    for (size_t g=0; g<n; g++) out[g]=rvec3_t(c[g], c[n+g], c[2*n+g]);
    return out;
}
template <class T, class Base> rvec3vec_t tSinglesGradientSampler<T,Base>::ProjectGradDM(const tDM_CD<T>& dm) const
{
    return Unpack(dm.Project(itsGrad).Vector());   // the density contracts its D into the gradient tensor
}
template <class T, class Base> rvec3vec_t tSinglesGradientSampler<T,Base>::ProjectGradField(const tChargeDensity<T>& cd) const
{
    return Unpack(itsGrad.ProjectField(cd).Vector());
}

template <class T, class Base> const rvec3vec_t& tSinglesGradientSampler<T,Base>::GradRhoPol(const tChargeDensity<T>* cd, const Spin& s) const
{
    assert(cd);
    if (s==Spin::None) throw std::logic_error("tSinglesGradientSampler::GradRhoPol: ask for a channel, not the total");
    if (cd->Version()!=itsGradVersion)
    {
        qchem::report::Timed timed("scf: XC-mesh grad-rho sampling (all iterations)");
        const tChargeDensity<T>* up=ChannelOf(cd, Spin::Up  );
        const tChargeDensity<T>* dn=ChannelOf(cd, Spin::Down);
        auto refuseMixed=[](const tChargeDensity<T>* c)
        {
            if (ExactSourceOf(c))
                throw std::logic_error("tSinglesGradientSampler: a density mixed as a field with its exact source retained (the cusp-deficit "
                                       "route) -- its correction is a G-series whose gradient no field gives yet, so a GGA on the point "
                                       "route needs rho[D] at the points (D-mixing); the policy layer should have chosen it");
        };
        if (auto* dmUp=dynamic_cast<const tDM_CD<T>*>(up))
        {
            auto* dmDn=dynamic_cast<const tDM_CD<T>*>(dn);
            if (!dmDn) throw std::logic_error("tSinglesGradientSampler: a D-backed Up channel beside a matrix-free Down channel");
            itsGradUp=ProjectGradDM(*dmUp);
            itsGradDn=ProjectGradDM(*dmDn);
        }
        else if (up && dn)
        {
            refuseMixed(up); refuseMixed(dn);
            itsGradUp=ProjectGradField(*up);
            itsGradDn=ProjectGradField(*dn);
        }
        else
        {   // spin-agnostic: grad rho_up = grad rho_down = grad rho / 2 (the HalfDensity rule RhoPol applies)
            if (auto dm=dynamic_cast<const tDM_CD<T>*>(cd)) itsGradUp=ProjectGradDM(*dm);
            else {refuseMixed(cd); itsGradUp=ProjectGradField(*cd);}
            for (auto& d : itsGradUp) d*=0.5;
            itsGradDn=itsGradUp;
        }
        itsGradVersion=cd->Version();
    }
    return s==Spin::Up ? itsGradUp : itsGradDn;
}

template <class T, class Base> const rvec_t& tSinglesGradientSampler<T,Base>::Rho(const tChargeDensity<T>* cd) const
{
    assert(cd);
    if (cd->Version()!=itsRhoSumVersion)
    {
        itsRhoSum=this->RhoPol(cd, Spin::Up)+this->RhoPol(cd, Spin::Down);
        itsRhoSumVersion=cd->Version();
    }
    return itsRhoSum;
}

template <class T, class Base> void tSinglesGradientSampler<T,Base>::WarmForDensity(const tChargeDensity<T>* cd, bool) const
{
    this->RhoPol(cd, Spin::Up); this->RhoPol(cd, Spin::Down);
    GradRhoPol(cd, Spin::Up);   GradRhoPol(cd, Spin::Down);
}

template <class T, class Base> template <class U>
hmat_t<U> tSinglesGradientSampler<T,Base>::MatrixVW(const tobs_t<U>* bs, const rvec_t& v, const rvec3vec_t& w) const
{
    const size_t n=this->NumPoints();
    if (v.size()!=n || w.size()!=n) throw std::invalid_argument("tSinglesGradientSampler::Matrix: v and w must be one entry per point");
    rvec_t c(3*n);
    for (size_t g=0; g<n; g++) {c[g]=w[g].x; c[n+g]=w[g].y; c[2*n+g]=w[g].z;}
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,T>&>(*bs);   // the DFT tier, as the LDA Matrix asks it
    hmat_t<U> H=this->template MatrixT<U>(bs, v);                                // v chi_i chi_j, through the LDA path
    H+=itsGrad.template Integrator<U>(orb).Adjoint(c);                           // + w . grad(chi_i chi_j), through the gradient tensor
    return H;
}
template <class T, class Base> hmat_t<T> tSinglesGradientSampler<T,Base>::Matrix(const tobs_t<T>* bs, const rvec_t& v, const rvec3vec_t& w) const
{
    return MatrixVW<T>(bs, v, w);
}

template class DeltaGradientProjector<double>;
template class DeltaGradientProjector<dcmplx>;
template class tSinglesGradientSampler<double>;                          // stage 2: the molecular point route
template class tSinglesGradientSampler<dcmplx, SinglesDensitySampler>;   // stage 3: the periodic point route
// The member templates the interface's inline faces call (a class instantiation does not reach them): each lineage's block scalars.
template const ScreenedMatrixIntegrator<double>& DeltaGradientProjector<double>::Integrator<double>(const BasisSet::Orbital_DFT_IBS<double,double>&) const;
template const ScreenedMatrixIntegrator<double>& DeltaGradientProjector<dcmplx>::Integrator<double>(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&) const;
template const ScreenedMatrixIntegrator<dcmplx>& DeltaGradientProjector<dcmplx>::Integrator<dcmplx>(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&) const;
template hmat_t<double> tSinglesGradientSampler<dcmplx, SinglesDensitySampler>::MatrixVW<double>(const tobs_t<double>*, const rvec_t&, const rvec3vec_t&) const;

} //namespace

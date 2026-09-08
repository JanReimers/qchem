// File: Hamiltonian/Internal/Imp/XCQuadrature_Singles.C  the SINGLES strategy: rho and H_xc both contracted through a cached Phi table.
//
// One implementation unit of module qchem.Hamiltonian.Internal.PWTerms.  Split 2026-09-08 out of a
// single 1213-line Imp/PWTerms.C (user: "PWTerms.C is huge, again doing too many things") into the
// interface-plus-many-Imp-units shape Internal/Terms.C has always had.  Helpers shared by more than
// one unit (NarrowExact, SampledField) live in the module INTERFACE's non-exported section, which is
// exactly what module-internal linkage is for -- they are visible to every unit of this module and to
// nothing outside it.
module;
#include <algorithm>   // std::min (the threaded quadrature's output-column blocking)
#include <cassert>
#include <complex>
#include <cstdlib>
#include <exception>   // std::exception_ptr (throw containment across the threaded Phi build)
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <optional>    // the conditionally-charged sub-buckets of the H_xc quadrature
#include <stdexcept>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.RunPolicy;   // theRunPolicy().XCFromDM() -- the declared XC-feed deviation (N5)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.G_FieldEvaluator;    // G_RasterTransform: the fit basis's FFT pair (RhoOnGrid, the BALL route)
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Parallel;                         // WorkerThreads (GPW_OMP_THREADS -- the XC-mesh table + quadrature loops)


namespace qchem::Hamiltonian
{

// ---- XC_SinglesQuadrature: the pair-shared mesh + Phi tables + per-serial rho ----------------------------------

XC_SinglesQuadrature::XC_SinglesQuadrature(fit_t fit, BasisSet::FitQuadrature quad)
    : itsFit(std::move(fit))
    , itsQuad(std::move(quad))
{
    // The fitter comes from the SAME Factory the molecular XC term uses; it inspects the basis and hands
    // back the delta fitter (R1.0 conformance).  A named lvalue because Factory takes its argument by
    // non-const reference, as its siblings do.
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> fb=itsFit;
    itsScalarFitter=Fitting::Factory(fb);
    assert(itsFit && "XC_SinglesQuadrature: the delta fit basis IS the quadrature -- it cannot be null");
    // The bundle's own invariants (fold partitions the mesh, sigmas need a fold, flags cover it) are
    // checked where the bundle becomes an object -- in the basis's ctor -- and the star-average announces
    // itself there too (providers self-report).  What IS this object's business is that the bundle it was
    // handed is the same one the basis holds: same length, so the fold's orbit indices and my coefficient
    // vectors index the same functions.
    assert((!itsQuad.mesh || itsQuad.mesh->size()==itsFit->GetNumFunctions()) &&
           "XC_SinglesQuadrature: the injected quadrature and the delta fit basis must be the SAME object "
           "-- one mesh point per fit function");
}

// THE FITTER'S PROJECTION FACE.  Since 2026-08-24 rho does NOT go to the fit basis: a projection is
// qcFitting's operation, and this strategy holds a fitter that already owns the fit basis's Phi handles
// for the adjoint direction (user).  So both directions of the SAME tensor now come off ONE object, which
// is the invariant this class's own header argues for one level up.  The cross-cast is the sanctioned
// "I want more" ask, made once per call on a face the delta fitter always carries.
const Fitting::ScalarProjector& XC_SinglesQuadrature::Projector() const
{
    auto* sp=dynamic_cast<const Fitting::ScalarProjector*>(itsScalarFitter.get());
    assert(sp && "XC_SinglesQuadrature: the delta scalar fitter must carry the projection face");
    return *sp;
}

// THE STAR-AVERAGE, applied to a coefficient vector over the delta basis (§6a W1).  It lives here, not on
// the fit basis, since 2026-08-24 (user): the operation is not a fitting question and the basis contributed
// only the fold -- which now arrives with the mesh, through the same factory out-parameter.  Both bodies are
// the free Symmetry::Lattice_3D algorithms, called directly.
void XC_SinglesQuadrature::Symmetrize(rvec_t& f) const
{
    if (itsQuad.fold.owner.empty()) return;              // free run: the projector is the identity
    Symmetry::Lattice_3D::SymmetrizeValues(itsQuad.fold, f);
}

// The MAGNETIC pair: with sigma tags the (rho,m) pair does NOT separate -- a Flip op maps rho_up onto
// rho_down, not onto itself -- so rho takes the plain orbit mean while m takes the chi-signed one, with m
// zeroed first at every point some Flip op fixes (where the exact projector annihilates it).  Without tags
// this is each channel on its own, which is bit-identical to what the removed base-class default did.
void XC_SinglesQuadrature::SymmetrizeSpin(rvec_t& rho, rvec_t& m) const
{
    if (itsQuad.fold.owner.empty() || itsQuad.sigmas.empty()) {Symmetrize(rho); Symmetrize(m); return;}
    for (size_t g=0; g<itsQuad.flipFixed.size(); ++g) if (itsQuad.flipFixed[g]) m[g]=0.0;
    Symmetry::Lattice_3D::SymmetrizeValues      (itsQuad.fold, rho);
    Symmetry::Lattice_3D::SymmetrizeValuesSigned(itsQuad.fold, itsQuad.sigmas, m);
}

// The quadrature questions go THROUGH the basis, and since 2026-08-23 in its own FUNCTION vocabulary:
// f is an expansion over the delta basis (c_g = f(r_g), which is what makes a pointwise functional
// applicable to it at all), so its integral is the coefficients dotted with the functions' own integrals,
// Integral(sum_g c_g delta_g) = sum_g c_g <delta_g|1> = sum_g c_g w_g.  Algebraically the old
// Integrate(values) and, written as this loop, the SAME summation order as qcMesh::Integrate was -- which
// is what keeps the pinned energies bit-unmoved.
const rvec_t& XC_SinglesQuadrature::FunctionIntegrals() const
{
    if (itsIntegrals.size()==0)
    {
        const vec_t<dcmplx> I=itsFit->Charge();
        itsIntegrals.resize(I.size());
        for (size_t a=0; a<I.size(); a++) itsIntegrals[a]=std::real(I[a]);
    }
    return itsIntegrals;
}
double XC_SinglesQuadrature::Integrate(const rvec_t& f) const
{
    const rvec_t& I=FunctionIntegrals();
    assert(f.size()==I.size() && "XC_SinglesQuadrature::Integrate: one coefficient per fit function");
    double s=0.0;
    for (size_t a=0; a<f.size(); a++) s+=I[a]*f[a];
    return s;
}
size_t XC_SinglesQuadrature::NumPoints() const {return itsFit->GetNumFunctions();}

// THE XC DM-rho REPAIR (doc/OpenWork.md, the factored-rho section).  Under rho-tilde mixing the density
// driving the Fock build is a G-space FIELD, so XC has to inverse-transform a truncated series at every
// mesh point -- 5.19 s per sampling on MnO (51% of the run) AND 8.5-18% of the points come back with
// rho<0, where the functionals guard `if (rho>0)` and silently contribute nothing.  But the mixer BUILT
// that field from a density matrix and now retains it (cDM_Sourced_CD), so the exact density is one
// cross-cast away: 0.042 s and rho>=0 by construction.  Hartree keeps the preconditioned field -- Poisson
// is LINEAR and diagonal in G -- while XC, a NONLINEAR POINTWISE functional, gets the cusp.  At the fixed
// point they agree, so this changes the SCF TRAJECTORY, not the answer.
// GPW_XC_DM_SOURCE=1 to arm it.  OPT-IN: a trajectory change must earn its place against banked recipes.
bool UseDMSource() {return theRunPolicy().XCFromDM();}
//! The exact density behind a field-backed one + THE DAMPING THAT FIELD APPLIED -- always together, because
//! taking the first without the second is the half-damped map (see cDM_Sourced_CD::EffectiveAlpha).
struct ExactSource
{
    const cDM_CD* cd=nullptr; double alpha=0.0;
    std::shared_ptr<const cChargeDensity> corr;   //!< N4: rho_mix - rho[D] as a field (null = wholesale route)
    explicit operator bool() const {return cd;}
};
// THE CUSP-DEFICIT ROUTE (doc/OpenWork.md N4):
//     rho_XC = rho[D]_exact + IFT[rho_mix - rho[D]]
// instead of the WHOLESALE replacement (damped rho[D]) GPW_XC_DM_SOURCE takes.  The difference is not
// cosmetic: this form's BAND-LIMITED content is identically rho_mix -- the array Hartree uses -- so Kerker's
// f(G) reaches XC unmodified and there is NO alpha_eff to choose.  A flat alpha_eff cannot reproduce that
// SHAPE, and measurement (2026-08-25) says destroying the selectivity costs the MnO magnetic basin
// (-61.403 -> -45.529, Eee 13.5 -> 29.0 Ha).
// NOT A FLAG HERE, BY DESIGN (user, 2026-08-25): whether the correction exists is decided when the MIXER is
// BUILT (SCFParams::XCCuspDeficit -> MakePeriodicMixer), so CP2K parity is a property of the constructed
// object and the plain Kerker mixer stays bit-identical.  This code just uses what it was handed: a
// non-null XCCorrection() means the run asked for the correction route.
ExactSource ExactSourceOf(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd)
{
    if (!cd) return {};
    auto* src=dynamic_cast<const qchem::ChargeDensity::cDM_Sourced_CD*>(cd);
    if (!src) return {};
    auto corr=src->XCCorrection();
    if (!corr && !UseDMSource()) return {};      // neither route armed
    return {src->DMSource().get(), src->EffectiveAlpha(), std::move(corr)};
}
// GPW_XC_DM_MIX overrides alpha_eff for CONTROLS only (=1 reproduces the undamped route); unset = use the
// mix's own.  GPW_XC_DM_BOOST scales it: alpha_eff came out ~0.20 on NaF and ~0.35 on MnO -- measured, but
// low against fractions those cells tolerate, and MnO converged 53 -> 39 iterations at boost 2.  NB f_K<=1
// forces alpha_eff<=alpha, so a boost >1/mean(f) puts XC ABOVE the mixer's own alpha -- defensible (XC's
// response kernel is finite at G->0, unlike Hartree's 4pi/G^2) but outside the preconditioner's bracket,
// hence a knob and not a default.
double DMSourceMixOverride()
{
    static const double a=[]{ const char* e=std::getenv("GPW_XC_DM_MIX"); return e ? std::atof(e) : -1.0; }();
    return a;
}
double DMSourceMixBoost()
{
    static const double b=[]{ const char* e=std::getenv("GPW_XC_DM_BOOST"); return e ? std::atof(e) : 1.0; }();
    return b;
}
//! Blend \a fresh into \a running at the mix's own alpha (outside (0,1) => passthrough, which is the right
//! bootstrap AND the right answer for an unmixed or overshooting step).  Non-negativity survives: a convex
//! combination of non-negative rasters is non-negative.
void DampXCChannel(rvec_t& running, const rvec_t& fresh, double alphaEff)
{
    const double ov=DMSourceMixOverride();
    const double a =(ov>=0.0) ? ov : DMSourceMixBoost()*alphaEff;
    static const bool trace=std::getenv("GPW_XC_ALPHA")!=nullptr;
    if (trace) std::cout<<"[XC alpha] alpha_eff="<<alphaEff<<(ov>=0.0?"  (OVERRIDDEN by GPW_XC_DM_MIX)":"")
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
//! Does \a cd carry a retained-D source (GPW_XC_DM_SOURCE / the N4 cusp deficit)?  A bool-returning
//! wrapper so a caller declared ABOVE ExactSource's definition can still ask.
bool HasExactSource(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd) {return bool(ExactSourceOf(cd));}

void ReportNegativeRho(const XC_Quadrature& q, const rvec_t& rho, const char* route)
{
    static const bool on=std::getenv("GPW_RHO_NEGATIVE")!=nullptr;
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
const rvec_t& XC_SinglesQuadrature::Rho(const cChargeDensity* cd) const
{
    assert(cd);
    // R2.9(i): the scalar and spin-resolved caches do not cross-invalidate (see the \warning on the
    // members).  An engine belongs to ONE xc/correlation pair, and a pair is either polarized or not, so
    // only one of the two routes is ever driven.  Pin it here rather than trusting the comment.
    assert(itsPolVersion==size_t(-1) && "XC_SinglesQuadrature: this engine already served RhoPol -- the scalar "
           "and spin-resolved rho caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsRhoVersion) return itsRho;
    itsRhoVersion=cd->Version();
    qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
    if (auto dm=dynamic_cast<const cDM_CD*>(cd))
    {
        itsRho=dm->ProjectOnto(Projector());   // the density contracts its D into the fitter's handles
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
            itsRho=ex.cd->ProjectOnto(Projector());
            const rvec_t c=Projector().Project(*ex.corr);
            assert(c.size()==itsRho.size() && "XC cusp-deficit: correction and rho must share the mesh");
            for (size_t g=0; g<itsRho.size(); ++g) itsRho[g]+=c[g];
            ReportNegativeRho(*this, itsRho, "cusp-deficit");
        }
        else
        {
            DampXCChannel(itsXCMix, ex.cd->ProjectOnto(Projector()), ex.alpha);
            itsRho=itsXCMix;   // the running mix lives in its OWN buffer -- see the \warning on itsXCMix
            ReportNegativeRho(*this, itsRho, "DM-source");
        }
    }
    else
    {
        itsRho=Projector().Project(*cd);   // non-DM (mixed rho-tilde / seed): the fitter projects the FIELD
        ReportNegativeRho(*this, itsRho, "matrix-free");
    }
    Symmetrize(itsRho);   // §6a W1: the injected quadrature's orbit-mean projector (no-op on a free run)
    return itsRho;
}

// The real-block ensure siblings (3c-3): build the real block's OWN typed table first (PhiR), then the
// shared sampling path -- the exact mirror of the complex ensureBlock argument.

// The spin-resolved sibling of Rho: the {up,down} PAIR is cached under ONE density serial (a polarized
// density's Version() forwards to its Up child -- a single scalar cache would hand the Up raster to the
// Down channel).  A cPolarized_CD answers per channel (each channel composite GEMMs its own D against the
// SHARED Phi tables); a spin-agnostic density (the seed) collapses to rho/2 per channel, so the first
// iterations run the exact unpolarized collapse (v^sigma(rho/2,rho/2)=v^P(rho)).
const rvec_t& XC_SinglesQuadrature::RhoPol(const cChargeDensity* cd, const Spin& s) const
{
    assert(cd);
    assert(itsRhoVersion==size_t(-1) && "XC_SinglesQuadrature: this engine already served the scalar Rho -- the "
           "two rho caches have no cross-invalidation, so one of them would go stale");
    assert(s!=Spin::None && "XC_SinglesQuadrature::RhoPol: ask for a channel, not the total");
    if (cd->Version()!=itsPolVersion)
    {
        itsPolVersion=cd->Version();
        qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
        if (auto pol=dynamic_cast<const ChargeDensity::cPolarized_CD*>(cd))
        {
            itsRhoUp=pol->GetChargeDensity(Spin::Up  )->ProjectOnto(Projector());
            itsRhoDn=pol->GetChargeDensity(Spin::Down)->ProjectOnto(Projector());
            ReportNegativeRho(*this, itsRhoUp, "DM(up)");
            ReportNegativeRho(*this, itsRhoDn, "DM(dn)");
        }
        else if (auto sr=dynamic_cast<const ChargeDensity::cSpinResolved_CD*>(cd))
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
            const ExactSource exUp=ExactSourceOf(sr->GetChannel(Spin::Up  ));
            const ExactSource exDn=ExactSourceOf(sr->GetChannel(Spin::Down));
            if (exUp && exDn)
            {
                // STALENESS GUARD (see itsSrcVersion): the field's serial just advanced, so the SOURCE's
                // must have too, or XC gets last iteration's D under a cache that believes it is fresh.
                if (exUp.cd->Version()==itsSrcVersion)
                    std::cerr<<"[XC DM-source] ** STALE: the mixed field advanced to serial "<<itsPolVersion
                             <<" but its retained density matrix did NOT (still "<<itsSrcVersion
                             <<") -- V_xc is being built from the PREVIOUS iteration's D."<<std::endl;
                itsSrcVersion=exUp.cd->Version();
                rvec_t up=exUp.cd->ProjectOnto(Projector());
                rvec_t dn=exDn.cd->ProjectOnto(Projector());
                if (exUp.corr && exDn.corr)
                {   // N4, per channel: each channel's own rho_mix - rho[D] (the polarized split gives the
                    // mixer one FourierMixCD per spin, so the corrections are already channel-private).
                    const rvec_t cu=Projector().Project(*exUp.corr), cd_=Projector().Project(*exDn.corr);
                    assert(cu.size()==up.size() && cd_.size()==dn.size());
                    for (size_t g=0; g<up.size(); ++g) {up[g]+=cu[g]; dn[g]+=cd_[g];}
                    itsRhoUp=up; itsRhoDn=dn;
                    ReportNegativeRho(*this, itsRhoUp, "cusp-deficit(up)");
                    ReportNegativeRho(*this, itsRhoDn, "cusp-deficit(dn)");
                }
                else
                {
                    DampXCChannel(itsXCMixUp, up, exUp.alpha);   // match the damping Hartree gets, so the map
                    DampXCChannel(itsXCMixDn, dn, exDn.alpha);   //   is not half-damped
                    itsRhoUp=itsXCMixUp; itsRhoDn=itsXCMixDn;
                    ReportNegativeRho(*this, itsRhoUp, "DM-source(up)");
                    ReportNegativeRho(*this, itsRhoDn, "DM-source(dn)");
                }
            }
            else
            {
            qchem::report::Timed seed("scf: XC-mesh rho sampling (matrix-free density)");
            itsRhoUp=Projector().Project(*sr->GetChannel(Spin::Up  ));
            itsRhoDn=Projector().Project(*sr->GetChannel(Spin::Down));
            ReportNegativeRho(*this, itsRhoUp, "matrix-free(up)");
            ReportNegativeRho(*this, itsRhoDn, "matrix-free(dn)");
            }
        }
        else
        {   // spin-agnostic seed: rho_up=rho_down=rho/2 (the molecular HalfDensity rule, cd85d13c)
            if (auto dm=dynamic_cast<const cDM_CD*>(cd))
                itsRhoUp=dm->ProjectOnto(Projector());
            else
                itsRhoUp=Projector().Project(*cd);
            itsRhoUp*=0.5;
            itsRhoDn=itsRhoUp;
        }
        {
            // Project the (ρ,m) PAIR over the injected quadrature: magnetic (σ tags -> ρ even, m odd with
            // the flip-fixed zeros), grey (each argument averaged independently) or free (a no-op).  The
            // three-way decision is ONE private method, not a branch spelled out at the call site.
            rvec_t rho = itsRhoUp + itsRhoDn;
            rvec_t m   = itsRhoUp - itsRhoDn;
            SymmetrizeSpin(rho, m);
            itsRhoUp = 0.5*(rho + m);
            itsRhoDn = 0.5*(rho - m);
        }
        // THE OBSERVABLE, reported where it is free (doc/OpenWork.md Step 0a).  This is the ONE place that
        // knows a NEW density has just been sampled on an atom-partitioned mesh, so the per-site integrated
        // moments cost a block sum over data already in hand -- and reporting them here means EVERY
        // polarized atom-centred run gets them, not just the one test that used to fake it with a point
        // probe.  Units: electrons (x mu_B for the magnetic moment).  Named partition, because until
        // Bader's zero-flux basins land the number is partition-dependent.
        EmitSiteMoments();
    }
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

// One line + one report entry per NEW density, from inside the serial-advance branch above.  Silent when
// the mesh carries no site partition (a uniform grid has no atomic basins to integrate over).
void XC_SinglesQuadrature::EmitSiteMoments() const
{
    const rvec_t mu=PartitionedMoments(rvec_t(itsRhoUp-itsRhoDn));
    if (mu.size()==0)
    {   // No site partition on this mesh -- legitimate for a uniform grid, a DEFECT for an atom-centred
        // one, and the difference used to be invisible: the instrument just printed nothing (it did so on
        // EVERY imposed run for as long as the invariant-mesh filter dropped the blocks).  Say which it is,
        // once, whenever the user asked for the moments.
        static bool said=false;
        if (std::getenv("QCHEM_SITE_MOMENTS") && !said)
        {
            said=true;
            std::cout<<"[site moments] UNAVAILABLE: the XC quadrature mesh carries no site blocks"
                     <<(itsQuad.mesh ? "" : " (no mesh injected at all)")
                     <<" -- an integrated site moment needs an atom-centred (Becke) mesh, and an "
                       "atom-centred mesh that lost its blocks is a defect, not a configuration."<<std::endl;
        }
        return;
    }
    double net=0.0, absSum=0.0;
    for (size_t a=0;a<mu.size();a++) { net+=mu[a]; absSum+=std::fabs(mu[a]); }
    if (absSum < 1e-8) return;                    // an unpolarized density has nothing to say
    qchem::report::json j;
    j["partition"]="Becke";                        // NOT canonical -- see the header's Bader note
    j["units"]="electrons";
    // Index loop, not the iterator-pair ctor: std cannot see Blaze's exported iterator op==/op!=
    // across the module boundary (CLAUDE.md "Includes & types").
    std::vector<double> v(mu.size());
    for (size_t a=0;a<mu.size();a++) v[a]=mu[a];
    j["mu"]=v;  j["net"]=net;
    qchem::report::EmitAt("scf", "siteMoments", j);
    if (std::getenv("QCHEM_SITE_MOMENTS"))
    {
        std::cout<<"[site moments] Becke-partitioned Integral w_A (rho_up-rho_dn) d3r [e]:";
        for (size_t a=0;a<mu.size();a++) std::cout<<"  "<<a<<":"<<mu[a];
        std::cout<<"   net="<<net<<std::endl;
    }
}

// The per-site INTEGRATED moment (see the header): mu_A = Integral w_A(r) [rho_up - rho_dn] d3r, in
// electrons.  Free -- RhoPol has already sampled both channels for this density serial (and cached them),
// and the mesh's weights already carry each site's Becke partition w_A, so this is a block sum over data
// in hand.  An unpolarized density gives exactly zero (rho_up == rho_dn by the HalfDensity collapse),
// which is the honest answer, not a special case.
rvec_t XC_SinglesQuadrature::SiteMoments(const cChargeDensity* cd) const
{
    assert(cd);
    const rvec_t& up=RhoPol(cd, Spin::Up);
    const rvec_t& dn=RhoPol(cd, Spin::Down);
    return PartitionedMoments(rvec_t(up-dn));      // empty when no partition was injected
}

// The ONE place the injected partition is read: Integral w_A f over each site block.  Not a fit-basis
// question and no longer asked of one -- the mesh arrives from the factory that built it for the basis,
// so this and the basis index the SAME points in the SAME order (one object, two collaborators).
rvec_t XC_SinglesQuadrature::PartitionedMoments(const rvec_t& f) const
{
    if (!itsQuad.mesh || itsQuad.mesh->NSites()==0) return rvec_t();
    assert(f.size()==itsQuad.mesh->size() && "XC_SinglesQuadrature: the injected quadrature and the fit "
           "basis must be the same object -- one field value per mesh point");
    return qcMesh::SiteIntegrals(*itsQuad.mesh, f);
}

// <i|v|j>: ASK THE BASIS.  It owns the points, the weights and the Phi table, so the whole quadrature
// -- Phi^dag diag(w v) Phi -- is its operation; this strategy only decides WHICH v to hand it.  That is
// what closed the last weight/coordinate escape (doc/CleanupCandidates.md R1.0 increment 2).
// <i|v|j> THROUGH THE FITTER, exactly as the molecular XC term does it (the Liskov conformance,
// doc/CleanupCandidates.md R1.0): fit the sampled field, then contract against this block.  The fit is
// re-done only when v CHANGES -- the Fock build calls this once per block with the same v, and DoFit
// on a delta basis is a copy, but re-copying per block would still be per-block work for nothing.
template <class U> hmat_t<U> XC_SinglesQuadrature::MatrixT(const tobs_t<U>* bs, const rvec_t& v) const
{
    bool same = itsFittedV.size()==v.size();
    for (size_t g=0; same && g<v.size(); g++) same = (itsFittedV[g]==v[g]);
    if (!same)
    {
        itsScalarFitter->DoFit(SampledField(v, NumPoints()));
        itsFittedV=v;
    }
    // TWO scalars on the face, and they DIFFER here (3c-3): the block is U, but the delta fit basis this
    // strategy holds is the run's periodic one, so the fit axis is dcmplx for both.  <U> alone would name
    // <U,U>, which for a real TRIM block is a face no fitter in the tree declares.  And the DFT cross-cast
    // is the XC strategy's, for the reason given on FitContraction: this layer is DFT by construction.
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>&>(*bs);
    return dynamic_cast<const Fitting::FitContraction<U,dcmplx>&>(*itsScalarFitter).Overlap(orb);
}
chmat_t XC_SinglesQuadrature::Matrix(const cobs_t* bs, const rvec_t& v) const {return MatrixT<dcmplx>(bs,v);}
rsmat_t XC_SinglesQuadrature::Matrix(const robs_t* bs, const rvec_t& v) const {return MatrixT<double>(bs,v);}

} //namespace

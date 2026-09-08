// File: Hamiltonian/Internal/Imp/XCQuadrature_Pair.C  the PAIR strategy: rho collocated through the 3-centre tensor, H by its raw adjoint.
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

// The grid-charge diagnostic's process-wide toggle.  It lives with the PAIR route because that is the only
// route that collocates rho onto a raster and can therefore report what the raster lost.
// ⚠ Process-wide MUTABLE state in a library -- flagged in doc/CleanupCandidates.md R1.0e as something the
// run report should own instead (theRunPolicy() already carries every other run-scoped switch).
bool& ReportGridCharge() { static bool on = false; return on; }

// ---- XC_PairQuadrature: rho by collocation, H by the SAME tensor's raw adjoint -------------------------

// It carries the fit basis (quadrature + collocation key + Overlap3C key) and, for the BALL fallback only,
// an ortho scalar fitter over it.  No grid announcement here: the fit basis self-reports at ITS
// construction, role-labeled (user ruling 2026-08-16).
XC_PairQuadrature::XC_PairQuadrature(fbs_t fb)
    : itsFitBasis(std::move(fb))
    , itsScalarFitter(Fitting::Factory(itsFitBasis))   // the ortho (G-space) scalar fit -- BALL route only
{
    assert(itsFitBasis);
}
XC_PairQuadrature::~XC_PairQuadrature() = default;   // itsScalarFitter's abstract type is complete here


// rho(r) on the raster for cd -- recomputed only on a new density serial, so the XC pair's two terms and
// their two energies share ONE collocation per iteration (it used to be one per TERM: the same raw
// collocation ran twice an iteration because each term owned its own copy of this).
// Forward declarations: these two helpers are defined below, beside the singles route that shares them.
// Declared rather than moved so the diff stays about the pair route (and so ONE definition still serves
// both, which is the point).
bool HasExactSource(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd);
void ReportNegativeRho(const XC_Quadrature& q, const rvec_t& rho, const char* route);

// ONE density object -> rho(r) on the raster.  RAW when it can collocate, BALL otherwise; extracted
// 2026-08-28 so the scalar route and the two spin channels cannot make that decision differently.
rvec_t XC_PairQuadrature::SampleOne(const cChargeDensity* cd, bool& isRaw) const
{
    assert(cd);
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "XC_PairQuadrature requires a FourierDensity (periodic) charge density");
    rvec_t rho=fd->GetRhoOnGrid(*itsFitBasis);
    isRaw=(rho.size()!=0);
    if (!isRaw)
    {
        auto* ge=dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get());
        assert(ge && "XC_PairQuadrature: the BALL route needs the fit basis's raster transforms");
        rho=ge->RhoOnGrid(fd->GetFourierDensity(*itsFitBasis));
    }
    return rho;
}

// ROUTE STABILITY (R2.16), in one place for both shapes -- see the declaration.
void XC_PairQuadrature::LatchRoute(const cChargeDensity* cd, bool isRaw) const
{
    if (!dynamic_cast<const ChargeDensity::cDM_CD*>(cd)) return;   // the seed cannot answer RAW; exempt it
    if (!itsRouteLatched) { itsRouteLatched=true; itsLatchedRaw=isRaw; return; }
    if (itsLatchedRaw!=isRaw)
        throw std::runtime_error(
            std::string("XC_PairQuadrature: the XC route changed mid-SCF (")
            + (itsLatchedRaw?"RAW -> BALL":"BALL -> RAW")
            + ").  These minimise DIFFERENT functionals -- BALL's ball-projected rho is non-variational "
              "-- so the optimiser would be chasing a moving target.  The route is a property of the "
              "orbital basis and must not change once the SCF is running.");
}

// THE SPIN-NATIVE SIBLING (2026-08-28).  Structurally the mirror of Refresh, once per channel -- and it
// walks the SAME two density shapes the singles route walks: a D-backed cPolarized_CD, and a matrix-free
// cSpinResolved_CD (the polarized seed, and the rho-tilde-mixed density on every Kerker/Pulay iteration).
void XC_PairQuadrature::RefreshPol(const cChargeDensity* cd) const
{
    assert(cd);
    assert(itsRhoVersion==size_t(-1) && "XC_PairQuadrature: this engine already served the scalar Rho -- the "
           "two rho caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsPolVersion) return;
    itsPolVersion=cd->Version();
    const cChargeDensity* up=nullptr;
    const cChargeDensity* dn=nullptr;
    if (auto pol=dynamic_cast<const ChargeDensity::cPolarized_CD*>(cd))
    {   up=pol->GetChargeDensity(Spin::Up); dn=pol->GetChargeDensity(Spin::Down); }
    else if (auto sr=dynamic_cast<const ChargeDensity::cSpinResolved_CD*>(cd))
    {   up=sr->GetChannel(Spin::Up);        dn=sr->GetChannel(Spin::Down); }

    qchem::report::Timed timed("scf: XC raw rho sampling (per channel)");
    if (!up || !dn)
    {
        // THE SPIN-AGNOSTIC SEED: rho_up = rho_dn = rho/2, the same HalfDensity rule the singles route
        // applies (cd85d13c).  Iteration 0's density has no channels to sample -- that is inherent to
        // seeding, not a gap here -- and its energy is discarded anyway.
        bool raw=false;
        itsRhoUp=SampleOne(cd, raw);
        itsRhoUp*=0.5;
        itsRhoDn=itsRhoUp;
        LatchRoute(cd, raw);                      // a no-op here: the seed carries no D, hence no latch
        // ★ THE SEED TAKES THE RAW ADJOINT EVEN WHEN ITS rho CAME FROM THE BALL ROUND TRIP, and that is
        // the exemption the route latch already grants in so many words: "a matrix-free density has no D
        // to collocate ... that is inherent to seeding, not a design choice, and iteration 0's energy is
        // discarded anyway".  Without it a SEEDED POLARIZED run on a real TRIM block has no adjoint at all
        // (the ball fit has no real contraction face -- see Matrix), which would put a grid restriction
        // back on a polarization property, i.e. re-make the conflation this change removed.
        // ⚠ What it costs: at iteration 0 ONLY, H_xc is not the exact derivative of the ball E_xc.  From
        // the first density-matrix-backed density both sides are RAW and the pairing is exact again --
        // and the latch makes that a checked property, not a hope.
        itsRhoIsRaw=raw;
        ReportNegativeRho(*this, itsRhoUp, "half-density seed");
        return;
    }
    // ⚠ NOT WIRED HERE: the DM-source repair and the N4 cusp deficit (GPW_XC_DM_SOURCE, XCCuspDeficit),
    // which the singles route applies per channel.  Both are default OFF.  Loud rather than silent,
    // because silently sampling the MIXED field where the caller asked for the retained D would be a
    // physics change wearing a performance change's clothes.
    if (HasExactSource(up) || HasExactSource(dn))
        throw std::logic_error("XC_PairQuadrature: GPW_XC_DM_SOURCE / XCCuspDeficit are not wired on the "
            "collocation (pair) route -- its rho comes from applyRaw, not from a projector, so the "
            "retained-D repair needs its own design.  Use the delta/singles quadrature for those flags.");
    bool rawUp=false, rawDn=false;
    itsRhoUp=SampleOne(up, rawUp);
    itsRhoDn=SampleOne(dn, rawDn);
    assert(rawUp==rawDn && "the two channels of one density must take the same XC route");
    LatchRoute(cd, rawUp);
    itsRhoIsRaw=rawUp;
    ReportNegativeRho(*this, itsRhoUp, "raw(up)");
    ReportNegativeRho(*this, itsRhoDn, "raw(dn)");
}

void XC_PairQuadrature::Refresh(const cChargeDensity* cd) const
{
    assert(cd);
    assert(itsPolVersion==size_t(-1) && "XC_PairQuadrature: this engine already served RhoPol -- the scalar "
           "and spin-resolved rho caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsRhoVersion) return;
    itsRhoVersion=cd->Version();
    // RAW-FIRST (doc/GPWPlan 0.5(f2)): a collocation-backed density answers with rho_DM(r) directly --
    // pointwise >= 0 for an aufbau D, so the rho>0 guard never bites and the grid calibration can relax
    // (the C=8 driver was the BALL path's Gibbs lobes).  A plane-wave density (or a matrix-free seed)
    // answers EMPTY and takes the ball round trip -- bit-identical to the pre-f2 path.  ONE decision,
    // in SampleOne, shared with the spin-resolved sibling.
    bool raw=false;
    itsRho=SampleOne(cd, raw);
    // THE SEED EXEMPTION, stated once for both shapes (see RefreshPol): a density with no D IS the seed,
    // it can only answer BALL, and it takes the RAW adjoint because the ball fit has no real contraction
    // face.  Iteration 0's energy is discarded and the latch below never fires on it, so the pairing is
    // exact from the first matrix-backed density onward.
    itsRhoIsRaw = raw;
    LatchRoute(cd, raw);
    // DIAGNOSTIC (env GPW_XCROUTE): which V_xc route fires this iteration -- RAW (applyRawAdjoint, FD-exact/
    // variational; gate GPW.RawXCConsistencyFD) vs BALL (DoFit/Overlap, non-variational under BallOnly;
    // gate GPW.XCPotentialConsistencyFD).  Answers whether GDM is fighting the ball non-variationality.
    if (std::getenv("GPW_XCROUTE"))
    {
        double rmin=1e300, rmax=-1e300;
        for (double r : itsRho) { rmin=std::min(rmin,r); rmax=std::max(rmax,r); }
        std::cout << "[xc route] " << (itsRhoIsRaw ? "RAW (applyRawAdjoint, variational)"
                                                   : "BALL (DoFit/Overlap, NON-variational under BallOnly)")
                  << "  rho grid=[" << std::scientific << std::setprecision(2) << rmin << ", " << rmax << "]"
                  << " npts=" << itsRho.size() << std::defaultfloat << std::endl;
    }
    if (ReportGridCharge())
    {
        // Grid charge vs analytic charge: the electrons LOST to grid truncation (high-G aliasing of rho).
        // == CP2K's "Electronic density on regular grids: <int rho> <error>" -- a controlled cutoff metric.
        const double qGrid=Integrate(itsRho);   // integral rho_grid d3r
        const double qDM  =cd->GetTotalCharge();                // Tr(D S) (analytic, ~ N)
        std::cout << "[grid charge] integral rho_grid=" << std::fixed << std::setprecision(6) << qGrid
                  << "  Tr(DS)=" << qDM
                  << "  lost=" << std::scientific << std::setprecision(3) << (qGrid-qDM)
                  << std::defaultfloat << std::endl;
        // XC-collapse diagnostic (doc/GPWPlan §0e step 2): the collocated rho's min/max/negative content.
        // If rho rings locally-negative near the sharp F peaks, GetEpsXc(rho<=0)=0 silently drops that
        // eps_xc*rho -- the +7 Ha Exc collapse.  Separates lead (c) [Gibbs ringing + guard: rho_min very
        // negative, big negCharge] from lead (b) [genuinely lower peaks: rho_max small, rho_min ~ 0].
        // NB the eps_xc estimate needs a functional and this object has none (the functional lives in the
        // TERM, which is the right place for it), so the lost-Exc column is reported by magnitude only.
        double rmin=1e300, rmax=-1e300; size_t nneg=0;
        rvec_t negOnly(itsRho.size());
        for (size_t q=0;q<itsRho.size();q++)
        {
            const double r=itsRho[q];
            rmin=std::min(rmin,r); rmax=std::max(rmax,r);
            negOnly[q] = r<0.0 ? r : 0.0;                                  // negative-charge density
            if (r<0.0) ++nneg;
        }
        std::cout << "[xc grid] rho_min=" << std::scientific << std::setprecision(3) << rmin
                  << " rho_max=" << rmax << " neg-frac=" << double(nneg)/double(itsRho.size())
                  << " negCharge=" << Integrate(negOnly)
                  << std::defaultfloat << std::endl;
    }
}

// Both go to the fit basis's RASTER face.  Not to the fit face: these two questions are about VOXELS (the
// array this strategy passes around is one value per voxel), and a plane-wave fit basis's own count is its
// {G} ball, which is a smaller and different number -- the exact confusion the 2026-08-23 "one fit-basis
// interface" pass removed by taking NumPoints/Integrate off FIT_SF_ABS.  Integral stays the raster's own
// (sum f)*Omega/N summation order, which the periodic energies are pinned to at 10 digits.
const BasisSet::G_RasterTransform& XC_PairQuadrature::Raster() const
{
    auto* ge=dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get());
    assert(ge && "XC_PairQuadrature: the pair route's fit basis is raster-backed by construction "
                 "(MakeXCQuadrature selects this strategy on exactly that capability)");
    return *ge;
}
double XC_PairQuadrature::Integrate(const rvec_t& f) const {return Raster().Integral(f);}
size_t XC_PairQuadrature::NumPoints() const {return Raster().RasterSize();}

const rvec_t& XC_PairQuadrature::Rho(const cChargeDensity* cd) const {Refresh(cd); return itsRho;}

const rvec_t& XC_PairQuadrature::RhoPol(const cChargeDensity* cd, const Spin& s) const
{
    assert(s!=Spin::None && "XC_PairQuadrature::RhoPol: ask for a channel, not the total");
    RefreshPol(cd);
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

// <i|v|j>, THE EXACT ADJOINT of whichever route Refresh took (which is why both live on one object):
//   RAW  -- the raw adjoint of the SAME tensor whose applyRaw produced itsRho (box-truncation per level +
//           the analytic gather), so H_xc == dE_xc/dD of the ONE raw discrete functional to machine
//           precision (gate: GPW.RawXCConsistencyFD).  No ball fit anywhere.
//   BALL -- fit v_xc on the {G} ball and contract: the legacy, non-variational pairing.
// Two-axis face (V1.1): the tensor follows TFit==dcmplx for both block scalars; narrow at the end.
template <class U> hmat_t<U> XC_PairQuadrature::MatrixT(const tobs_t<U>* bs, const rvec_t& v) const
{
    assert((itsRhoVersion!=size_t(-1) || itsPolVersion!=size_t(-1))
           && "XC_PairQuadrature::Matrix before Rho/RhoPol: the adjoint must follow the collocation that "
              "fixed the route");
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>&>(*bs);   // genuine "is it?" cross-cast (throws)
    const Projector3<dcmplx>& g=orb.Overlap3C(*itsFitBasis);
    // ★ THE ADJOINT FOLLOWS THE LINEAGE'S CAPABILITY; rho follows the DENSITY's.  They differ on exactly
    // one iteration and the design already grants it: a matrix-free SEED has no D to collocate, so it can
    // only sample BALL, while the route latch exempts it and iteration 0's energy is discarded anyway.
    // From the first matrix-backed density a GPW lineage always answers RAW, so the pairing is exact from
    // there on -- and the latch makes that a CHECKED property rather than a hope.
    //
    // ⚠ WHY THE TEST IS THE CAPABILITY AND NOT THE SEED FLAG.  A plane-wave lineage never sets applyRaw at
    // all, so it takes the ball fit throughout -- as it always has.  A GPW lineage always has the raw
    // adjoint.  So `does this lineage have one?` decides it, and it can only be asked HERE, with `bs` in
    // hand.  Trying to decide it back in Refresh needs a "is this the seed?" flag, and the polarized seed
    // showed why that is the wrong question: PolarizedSeedCD IS spin-resolved, so it HAS channels and looks
    // matrix-backed from the parent, while its channels carry no D.
    //
    // ⚠ AND WITHOUT THIS a seeded POLARIZED run on a real TRIM block has no adjoint at all -- the ball fit
    // has no FitContraction<double,dcmplx> face -- which would put a grid restriction back on a
    // polarization property, i.e. re-make the conflation this change removed.
    if (g.applyRawAdjoint)
        return NarrowExact<U>(g.applyRawAdjoint(v));
    if constexpr (std::is_same_v<U,dcmplx>)
    {
        itsScalarFitter->DoFit(SampledField(v, NumPoints()));
        // The COMPLEX contraction face (ISP): this fitter's raster basis serves Bloch blocks on a complex
        // fit axis.  The DFT cross-cast is THIS strategy's to make (2026-08-24): the contraction face takes
        // the DFT block, and an XC strategy is DFT-specific by construction -- unlike the generic term
        // machinery it is reached through, which must stay method-neutral.
        return dynamic_cast<const Fitting::FitContraction<dcmplx,dcmplx>&>(*itsScalarFitter).Overlap(orb);
    }
    else
        // ⛔ A REAL TRIM BLOCK ON THE BALL FIT IS STILL UNWIRED, and now we know exactly why rather than
        // "nothing reaches it": the ortho scalar fitter implements FitContraction<dcmplx,dcmplx> and NOT
        // the <double,dcmplx> face, so the cast is a std::bad_cast (verified 2026-08-28).  The face is
        // expressible -- FitContraction is templated on the block scalar -- so this is a hole in the
        // FITTING layer, not a fact about XC.  The seed exemption below is what keeps it unreachable.
        throw std::logic_error("XC_PairQuadrature: a real TRIM block on the legacy ball-fit XC route is not "
                               "wired -- the ortho scalar fitter has no FitContraction<double,dcmplx> face. "
                               "Use the raw collocated feed (the default) or the delta quadrature");
}
chmat_t XC_PairQuadrature::Matrix(const cobs_t* bs, const rvec_t& v) const {return MatrixT<dcmplx>(bs,v);}
rsmat_t XC_PairQuadrature::Matrix(const robs_t* bs, const rvec_t& v) const {return MatrixT<double>(bs,v);}

// CAPABILITY DECIDES (doc/OpenWork.md): a delta basis carries points and nothing else -> singles; a
// raster-backed one carries the FFT transforms and keys the 3-centre tensor -> pair.  One decision, taken
// once, and latched for the run by the simple fact that the Hamiltonian builds this object once.
std::shared_ptr<const XC_Quadrature>
MakeXCQuadrature(const std::shared_ptr<const BasisSet::cFIT_SF_ABS>& fb,
                 BasisSet::FitQuadrature quad)
{
    assert(fb);
    // CAPABILITY decides.  A basis that carries the {r}<->{G} transforms is raster-backed, so its
    // collocation pair (Overlap3C's applyRaw/applyRawAdjoint) exists and the PAIR route is available --
    // and preferred, being the production GPW path.  Anything else can only be contracted through a Phi
    // table: SINGLES.  Note this asks what the basis CAN do, never what it IS.
    if (dynamic_cast<const BasisSet::G_RasterTransform*>(fb.get()))
        return std::make_shared<const XC_PairQuadrature>(fb);
    return std::make_shared<const XC_SinglesQuadrature>(fb, std::move(quad));
}

} //namespace

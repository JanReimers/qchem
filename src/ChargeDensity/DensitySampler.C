// File: ChargeDensity/DensitySampler.C  The SCF iteration's rho on a fit basis's sampling axis --
// the ABSTRACT face and its factory.  Nothing else.
//
// ★★ WHY IT LIVES IN qcChargeDensity (R1.0e, settled 2026-09-10).  Its whole interface takes a
// `const cChargeDensity*` -- Rho, RhoPol, WarmForDensity -- so the DAG decides the question
// before taste gets a vote: qcChargeDensity sits ABOVE qcFitting, which means qcFitting could not host this
// even if it wanted to (it cannot import qchem.ChargeDensity).  And it does not want to: every reference to
// this engine in qcFitting is a COMMENT citing it as precedent.  qcChargeDensity, by contrast, already
// describes how its own classes serve this consumer in five places (IrrepCD::ProjectOnto, the RhoPol
// cross-cast, DensityMixer, FourierDensity).
//
// ⇒ It was never an XC-library citizen.  It touches a functional ZERO times; a Hartree term or a +U
// projector would want the same object, which is why it stopped being called `XC_Quadrature` (2026-09-09)
// and why it stopped living in qcHamiltonian (here).  Moved BEFORE DFT+U is written against it, because
// moving it afterwards would churn +U too.
//
// ★ THE SPLIT IS THE PROJECT'S STANDING PATTERN (user, 2026-09-10: *"abstract interface + factory and
// concrete imp internal is a consistent pattern for this project.  Give the clients what they need and
// nothing more."*):
//   - HERE, public:  the abstract `DensitySampler` + `MakeDensitySampler`.  That is the whole client surface.
//   - `qchem.ChargeDensity.Internal.DensitySampler`: the two concrete strategies.  A client never names one.
// The factory is DEFINED in this module's own Imp unit, which may import the Internal module -- so the
// concretes stay invisible to everyone above.
module;
#include <cassert>
#include <cstddef>
#include <complex>
#include <stdexcept>
#include <iosfwd>
#include <map>
#include <optional>   // GroupResidual (R2.1)
#include <memory>
#include <string>
#include <vector>   // SinglesDensitySampler sigmas/flipFixed (Shubnikov S3)
export module qchem.ChargeDensity.DensitySampler;
import qchem.BasisSet.Orbital_DFT_IBS;      // the fit-basis faces
import qchem.BasisSet.Fit_Types;            // FitQuadrature
import qchem.ChargeDensity.Types;           // tobs_t/cobs_t/robs_t -- this library has its OWN, identical to the
                                            // qcHamiltonian typedefs the engine used to borrow (R1.0e, 2026-09-10)
                                            // with no Hamiltonian dependency of its own; it moves with the engine
import qchem.ChargeDensity;
export import qchem.ChargeDensity.TransitionDensity;   // Sample: δρ of a linear response (R2)
import qchem.Mesh;                          // qcMesh::Mesh/MeshParams (the quadrature the engine integrates on)
export import qchem.Mesh.Integrator;        // qcMesh::MatrixAdjoint -- the ONE face this engine names
import qchem.Types;
import qchem.PointSampling;                  // SamplingTarget (FourierTarget's slot)

export namespace qchem::ChargeDensity
{

// The density names this engine consumes, pulled in EXPLICITLY rather than inherited from
// qchem.Hamiltonian.  The engine must not import the term face -- that is the dependency the extraction
// exists to break -- so it states for itself what it takes.  This list IS the engine's coupling to
// qcChargeDensity, and it is what makes qcFitting an illegal home (see the header).
using ChargeDensity::cChargeDensity;
using ChargeDensity::cDM_CD;
using ChargeDensity::tChargeDensity;

//! Process-wide diagnostic toggle (default OFF).  When true,
//! \c PairDensitySampler::Refresh emits a one-line report each time it (re)collocates the density: the grid-integrated
//! charge \f$\int\rho_{\text{grid}}\f$, the analytic charge \f$\mathrm{Tr}(DS)\f$, and their difference -- the
//! CHARGE LOST TO GRID TRUNCATION (== CP2K's "Electronic density on regular grids: <int> <error>" readout).
//! A cheap, controlled number for "is the density cutoff high enough" (see doc/GPWPlan.md \S0).  Flip in place:
//! `qchem::ChargeDensity::ReportGridCharge() = true;`.
bool& ReportGridCharge();

//! \brief WHICH \f$\rho\f$ a sampled raster came from -- the pin-28 audit (doc/Pins.md pin 28: XC sees \f$\rho[D]\f$,
//! never a fitted \f$\rho\f$).  One value per sampling; the XC term records the one its \f$V_{xc}\f$ build saw.
enum class XCRhoSource
{
    DM,           //!< \f$\rho[D]\f$ of a density matrix, pointwise
    CuspDeficit,  //!< N4: \f$\tilde\rho_{mix}+(\rho[D]-\tilde\rho_{fit})\f$ -- pin 28's exception 2, in its legal form
    DMSource,     //!< the wholesale retained-D route (policy.xcFromDM): a damped mix of \f$\rho[D]\f$
    Shadow,       //!< the raster route's RAW shadow of a mixed density: \f$\rho[D]\f$ collocated, Kerker-mixed on the raster
    //! An EXACT fit of \f$\rho[D]\f$: rho always goes through the fit (user's definition, 2026-10-26: that keeps the SCF loop
    //! basis-independent), and for a plane-wave orbital basis whose fit ball holds every pair difference (\c G_OrbitalPairBall,
    //! >= 4 Ecut) the fit is Parseval-exact -- \f$\rho[D]\f$ to roundoff.  Legal under pin 28.
    ExactFit,
    Fitted,       //!< an INEXACT fit / projected field (the seed, a bare \f$\tilde\rho_{mix}\f$, a truncated ball) -- legal only for the seed
    //! \f$P\rho\f$ BY STATED CHOICE (\c xc.rho "fit" on an inexact fit, V-XC-ONE-TERM step 3): the variational projected cell, where XC sees
    //! the fit's projection of \f$\rho[D]\f$ because the user asked for that functional.  An explicit, reported departure from pin 28 -- so
    //! it is NOT \c Fitted, which marks a fitted \f$\rho\f$ nobody asked for.
    Projected,
};
//! The short label the SCF trace's "XC ρ" column prints for \a s.
inline const char* ShortName(XCRhoSource s)
{
    switch (s)
    {
        case XCRhoSource::DM:          return "DM";
        case XCRhoSource::CuspDeficit: return "cusp";
        case XCRhoSource::DMSource:    return "DMsrc";
        case XCRhoSource::Shadow:      return "shadow";
        case XCRhoSource::ExactFit:    return "fit=";
        case XCRhoSource::Fitted:      return "fitted";
        case XCRhoSource::Projected:   return "Pρ";
    }
    assert(false && "ShortName: unhandled XCRhoSource");
    return "?";
}

//! \brief THE XC QUADRATURE: ONE OBJECT WITH TWO ADJOINT-PAIRED FACES -- \f$\rho\f$ at the points from
//! the density, and \f$\langle i|v|j\rangle\f$ from a field at those SAME points.
//!
//! WHY THEY ARE ONE OBJECT AND NOT TWO AXES (doc/OpenWork.md, "separation of concerns in the XC terms").
//! \f$H=\partial E/\partial D\f$ holds only if the integrate-back is the EXACT ADJOINT of the collocation
//! ON THE SAME TRUNCATED OPERATOR -- what the \c GPW.RawXCConsistencyFD gate enforces.  Split the two
//! across separate owners and a mismatch becomes EXPRESSIBLE: measured 2026-08-22, an unscreened singles
//! \f$\rho\f$ paired with a screened pair \f$H\f$ sent Si from 14 to 60 iterations and moved E by 35 μHa.
//! With both faces on one type, the pairing is a class invariant instead of a convention, and "GDM just
//! works" stops resting on discipline.
//!
//! The two implementations below are a COST/TRUNCATION STRATEGY, not a semantics: on the same points both
//! compute \f$H_{ij}=\sum_g w_g v(r_g)\chi_i(r_g)\chi_j(r_g)\f$, differing only in evaluation order and
//! hence in what each truncates (the pair route's ε-screening and multigrid boxes).  Neither always wins --
//! pair scales with the SCREENED pair count, singles with \f$n_{pts}n^2\f$ -- so which one runs is decided
//! by \c MakeDensitySampler from the fit basis's capabilities, once, and LATCHED for the run (switching
//! mid-SCF would change the truncated operator, i.e. the functional).
//! \brief The scalar-free half of the XC quadrature: what it integrates and over how many points.  Shared by both scalar
//! types (V-XC-ONE-TERM, 2026-10-31): a molecule's real density and a crystal's Bloch density sample onto the same kind of
//! point set, and nothing here knows which.
class DensitySamplerBase
{
public:
    virtual ~DensitySamplerBase() = default;
    //! \f$\int f\,d^3r\f$ for a field sampled at MY points -- the \f$E_{xc}\f$ quadrature.  A term hands
    //! back a value array and never learns where the points are (nor which kind of mesh they came from).
    //! \note POINT vocabulary is correct HERE and nowhere below it: this face IS a quadrature (its whole
    //! job is \f$\rho\f$ at points and the adjoint back), which is exactly what the fit BASIS is not.  How
    //! each strategy answers differs accordingly -- the \f$\delta\f$ one dots the coefficients with its
    //! functions' integrals (\c FIT_SF_ABS::Integrals), the raster one uses the raster's uniform rule.
    virtual double Integrate(const rvec_t& f) const=0;
    //! How many points I sample at -- for reporting only (a term's \c Write line).
    virtual size_t NumPoints() const=0;
    //! The source of the \f$\rho\f$ I have CACHED -- i.e. what the last \c Rho / \c RhoPol call served, hit or miss.
    //! Client: the XC term, right after it sampled for a \f$V_{xc}\f$ build, so the SCF trace can show what the
    //! Fock matrix's XC saw each iteration (pin 28).  Asked of the term's own build rather than logged here, because
    //! this engine also samples for the energy pass and the residual, which would bury the one sampling that matters.
    virtual XCRhoSource RhoSource() const=0;
    //! \brief R2.1 (doc/SCFRobustnessPlan.md P2): the density this sampler last cached, measured against the imposed group
    //! BEFORE its own projector -- \f$\max|f-Pf|/\max|\rho|\f$ over \f$\rho\f$ (and \f$m\f$ on a polarized run).  What XC reads is
    //! \f$Pf\f$; the number says how far from invariant the density it was HANDED was.  0 on a free run (P is the identity);
    //! nullopt from a sampler that does not project at all, which therefore cannot answer.
    virtual std::optional<double> GroupResidual() const=0;
    //! \brief \f$\int w_A f\,d^3r\f$ per SITE BLOCK of my quadrature, for a field sampled at MY points -- the
    //! atom-partitioned sibling of \c Integrate.  Default EMPTY: a quadrature with no atomic partition (a
    //! uniform raster) has no basins to integrate over, and the caller must ask rather than assume.
    //! A QUADRATURE question and nothing more (R1.0h): the OBSERVABLE built on it -- the integrated site
    //! moment \f$\mu_A\f$ of \f$\rho_\uparrow-\rho_\downarrow\f$ -- is the spin-native XC term's to
    //! compute and \c ChargeBreakdown's to carry; this engine no longer knows what it is integrating.
    virtual rvec_t SiteIntegrals(const rvec_t& f) const {return rvec_t();}
    //! \brief The sampling slot of MY points (the rho-tilde sampler, step 2): the accuracy the run stated for this quadrature and the ONE
    //! tolerance tau a Fourier series sampled here is held to.  Client: the XC term, which derives tau from the seed density (it alone
    //! has v_xc and f_xc in hand) and reports what was delivered.  Default null: a quadrature that samples no Fourier series by NUFFT
    //! (the raster pair route) has no slot to size.
    virtual SamplingTarget* FourierTarget() const {return nullptr;}
};

//! \brief The XC quadrature for densities and blocks of scalar type \a T -- the faces \c Vxc_Quadrature<T> drives.
//! \c double for the molecular / atom path, \c dcmplx for the periodic one; the periodic sampler adds the response and the
//! real-TRIM faces below (\c DensitySampler).  (V-XC-ONE-TERM: until 2026-10-31 this face was complex-only, which is why a
//! molecule could not run the one XC term.)
template <class T> class tDensitySampler
    : public virtual DensitySamplerBase
{
public:
    //! \brief \f$\rho(r_g)\f$ for \a cd's current serial, cached across the XC pair.
    //! (No "ensure this block is tabled first" hint any more: the density now asks the QUADRATURE for each
    //! of its own blocks' tables, so there is no first-pass gap for a caller to plug -- 2026-08-22.)
    virtual const rvec_t& Rho(const tChargeDensity<T>* cd) const=0;
    //! Spin channel \f$\rho_\sigma(r_g)\f$ -- the SPIN-NATIVE sibling of \c Rho (§4 tier 4b).  Not every
    //! quadrature can answer it (the pair route has no per-spin collocation): those THROW, and the
    //! Hamiltonian's Auto rule keeps a polarized run off them.
    virtual const rvec_t& RhoPol(const tChargeDensity<T>* cd, const Spin& s) const=0;
    //! \f$\langle i|v|j\rangle=\sum_g w_g\,\overline{\chi_i(r_g)}v_g\chi_j(r_g)\f$ -- the EXACT ADJOINT of
    //! whatever route \c Rho took, weights included (a caller passes the bare field \f$v\f$).
    virtual hmat_t<T> Matrix(const tobs_t<T>* bs, const rvec_t& v) const=0;
    //! THE DERIVATIVE BESIDE \c Matrix (doc/GeometryOptimizationPlan.md §13.2; lane B 2026-11-09): per atom of \a st,
    //! Re Tr(D ∂_R Matrix(bs, v)) at the FIXED field \a v -- the XC term's force at fixed D, since \f$\partial E_{xc}/\partial\rho\f$
    //! at the points IS v and this quadrature's \c Matrix is the exact adjoint of its \c Rho.  CLIENT: the XC term's force face.
    //! The default REFUSES by name: a strategy that has not stated its derivative has none (a moving mesh, a table on fixed points).
    //! \c CanPullbackMatrix answers the same question without the work (CLIENT: the XC term's MissingGradient, at setup).
    virtual bool CanPullbackMatrix(const tobs_t<T>*) const {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const tobs_t<T>* bs, const rvec_t&, const hmat_t<T>&, const Structure*) const
    {return Outcome<rvec3vec_t,AD::NoGradient>::Fail({"the XC quadrature's strategy gives no pullback of its matrix on the block "+bs->Name()});}
    //! \brief Pre-warm this engine's per-density caches for \a cd -- the EAGER REFRESH PHASE
    //! (doc/OpenWork.md item **KP**).  \a polarized picks which shape to warm, because the two are
    //! mutually exclusive on one engine (see the cross-invalidation warning on both implementations: an
    //! engine answers \c Rho or \c RhoPol for a run, never both).
    //!
    //! It is expressed as ONE call rather than "the term calls Rho/RhoPol itself" so that the WARMING and
    //! the SHAPE RULE stay in the engine that owns the caches; a term that reached in by name would have
    //! to know the exclusivity rule too.
    virtual void WarmForDensity(const tChargeDensity<T>* cd, bool polarized) const
    {
        if (polarized) { RhoPol(cd, Spin::Up); RhoPol(cd, Spin::Down); }
        else           { Rho(cd); }
    }
};

//! \brief The PERIODIC quadrature: \c tDensitySampler<dcmplx> plus the faces only a Bloch run has -- the linear-response
//! transition pair (R2/R3) and the real TRIM block (Step 3c).  Its name is the one every periodic client already spells.
class DensitySampler
    : public virtual tDensitySampler<dcmplx>
{
public:
    using tDensitySampler<dcmplx>::Matrix;
    //! \brief \f$\delta\rho_\sigma(r_g)\f$ of a TRANSITION density at my points (doc/LinearResponsePlan.md R2) --
    //! through the SAME route and projector as \c Rho, so a response samples on exactly the quadrature its
    //! ground state did.  \a s = \c Spin::None asks for the total (the folded doublet).
    //! UNCACHED, because it must not evict the SCF iterate's \f$\rho\f$ (the ground-state caches stay warm);
    //! and NEVER SYMMETRIZED, because a perturbation may break the imposed symmetry and its response must be
    //! free to (the fold's orbit mean is a property of the ground state, not of δρ).
    //! ★ COMPLEX, and in THIS QUADRATURE'S OWN REPRESENTATION of a Bloch-q function (R3 step 3): the point route
    //! returns δρ at its points, the raster route the periodic part \f$e^{-iq\cdot r}\delta\rho\f$ at its raster
    //! points.  Either way a PERIODIC pointwise kernel (\f$f_{xc}[\rho_0]\f$) may multiply it, and the pair
    //! \c Matrix below is its exact adjoint -- which is all a term needs to know.  Real at q = 0 for a Hermitian δD.
    virtual cvec_t Sample(const TransitionDensity<dcmplx>& delta, const Spin& s) const=0;
    //! \brief The transition ADJOINT of \c Sample: \f$h_{ij}=\langle\chi_i^{k+q}|\,\tilde v\,|\chi_j^k\rangle\f$ for a
    //! field \a v in \c Sample's representation, bra (k+q) x ket (k).  \a q is the wave vector the field was sampled
    //! at (the raster route gathers a periodic part and needs it; the point route does not).
    virtual mat_t<dcmplx> Matrix(const cobs_t* bra, const cobs_t* ket, const cvec_t& v, const rvec3_t& q) const=0;
    //! The REAL-BLOCK sibling (Step 3c): a real TRIM block's quadrature runs in REAL arithmetic.
    virtual rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const=0;
    using tDensitySampler<dcmplx>::MakeMatrixPullback;
    using tDensitySampler<dcmplx>::CanPullbackMatrix;
    virtual bool CanPullbackMatrix(const robs_t*) const {return false;}
    //! The REAL-BLOCK sibling of \c MakeMatrixPullback (the derivative beside the real \c Matrix).  Default: a named refusal.
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeMatrixPullback(const robs_t* bs, const rvec_t&, const hmat_t<double>&, const Structure*) const
    {return Outcome<rvec3vec_t,AD::NoGradient>::Fail({"the XC quadrature's strategy gives no pullback of its matrix on the real block "+bs->Name()});}
};

//! \brief Pick the assembly strategy for \a fb -- CAPABILITY decides, and the answer is fixed for the run.
//!
//! A δ fit basis carries points and nothing else, so it can only be contracted
//! through a Φ table: SINGLES.  A raster-backed one additionally carries the FFT transforms and keys the
//! orbital's 3-centre tensor, so the PAIR route is available and is chosen -- it is the production GPW
//! path and the one whose screening pays on large cells.  There is no extra input to supply: both routes
//! are already functions of (orbital basis, fit basis).
//! The fit basis a sampler is built over -- named HERE because it is the FACTORY's parameter type, and a
//! client that has to spell it should not have to reach into a concrete strategy for the spelling.
typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fitbasis_t;

std::shared_ptr<const DensitySampler>
MakeDensitySampler(const fitbasis_t& fb, BasisSet::FitQuadrature quad={});

//! \brief The REAL path's strategy (V-XC-ONE-TERM, 2026-10-31): a molecule's or an atom's δ basis -- points and nothing else,
//! so SINGLES.  A non-orthogonal (A1 Gaussian) fit is a projection: the PROJECTED strategy, \f$P\rho\f$ evaluated on \a quad's mesh
//! (step 3, 2026-11-10) -- the \c XCρ::Fit cell, which the caller states.
std::shared_ptr<const tDensitySampler<double>>
MakeDensitySampler(const std::shared_ptr<const BasisSet::rFIT_SF_ABS>& fb, BasisSet::FitQuadrature quad={});

} //namespace


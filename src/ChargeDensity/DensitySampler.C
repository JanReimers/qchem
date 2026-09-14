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
#include <memory>
#include <string>
#include <vector>   // SinglesDensitySampler sigmas/flipFixed (Shubnikov S3)
export module qchem.ChargeDensity.DensitySampler;
import qchem.BasisSet.Orbital_DFT_IBS;      // the fit-basis faces + FitQuadrature
import qchem.ChargeDensity.Types;           // tobs_t/cobs_t/robs_t -- this library has its OWN, identical to the
                                            // qcHamiltonian typedefs the engine used to borrow (R1.0e, 2026-09-10)
                                            // with no Hamiltonian dependency of its own; it moves with the engine
import qchem.ChargeDensity;
import qchem.Mesh;                          // qcMesh::Mesh/MeshParams (the quadrature the engine integrates on)
export import qchem.Mesh.Integrator;        // qcMesh::MatrixAdjoint -- the ONE face this engine names
import qchem.Blaze;                         // blazem::NarrowExact (the real-TRIM narrow, promoted to qcMath 2026-09-08)
import qchem.Types;

export namespace qchem::ChargeDensity
{

// The density names this engine consumes, pulled in EXPLICITLY rather than inherited from
// qchem.Hamiltonian.  The engine must not import the term face -- that is the dependency the extraction
// exists to break -- so it states for itself what it takes.  This list IS the engine's coupling to
// qcChargeDensity, and it is what makes qcFitting an illegal home (see the header).
using ChargeDensity::cChargeDensity;
using ChargeDensity::cDM_CD;

//! Process-wide diagnostic toggle (default OFF).  When true,
//! \c PairDensitySampler::Refresh emits a one-line report each time it (re)collocates the density: the grid-integrated
//! charge \f$\int\rho_{\text{grid}}\f$, the analytic charge \f$\mathrm{Tr}(DS)\f$, and their difference -- the
//! CHARGE LOST TO GRID TRUNCATION (== CP2K's "Electronic density on regular grids: <int> <error>" readout).
//! A cheap, controlled number for "is the density cutoff high enough" (see doc/GPWPlan.md \S0).  Flip in place:
//! `qchem::ChargeDensity::ReportGridCharge() = true;`.
bool& ReportGridCharge();

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
class DensitySampler
{
public:
    virtual ~DensitySampler() = default;
    //! \f$\int f\,d^3r\f$ for a field sampled at MY points -- the \f$E_{xc}\f$ quadrature.  A term hands
    //! back a value array and never learns where the points are (nor which kind of mesh they came from).
    //! \note POINT vocabulary is correct HERE and nowhere below it: this face IS a quadrature (its whole
    //! job is \f$\rho\f$ at points and the adjoint back), which is exactly what the fit BASIS is not.  How
    //! each strategy answers differs accordingly -- the \f$\delta\f$ one dots the coefficients with its
    //! functions' integrals (\c FIT_SF_ABS::Integrals), the raster one uses the raster's uniform rule.
    virtual double Integrate(const rvec_t& f) const=0;
    //! How many points I sample at -- for reporting only (a term's \c Write line).
    virtual size_t NumPoints() const=0;
    //! \brief \f$\rho(r_g)\f$ for \a cd's current serial, cached across the XC pair.
    //! (No "ensure this block is tabled first" hint any more: the density now asks the QUADRATURE for each
    //! of its own blocks' tables, so there is no first-pass gap for a caller to plug -- 2026-08-22.)
    virtual const rvec_t& Rho(const cChargeDensity* cd) const=0;
    //! Spin channel \f$\rho_\sigma(r_g)\f$ -- the SPIN-NATIVE sibling of \c Rho (§4 tier 4b).  Not every
    //! quadrature can answer it (the pair route has no per-spin collocation): those THROW, and the
    //! Hamiltonian's Auto rule keeps a polarized run off them.
    virtual const rvec_t& RhoPol(const cChargeDensity* cd, const Spin& s) const=0;
    //! \f$\langle i|v|j\rangle=\sum_g w_g\,\overline{\chi_i(r_g)}v_g\chi_j(r_g)\f$ -- the EXACT ADJOINT of
    //! whatever route \c Rho took, weights included (a caller passes the bare field \f$v\f$).
    virtual chmat_t Matrix(const cobs_t* bs, const rvec_t& v) const=0;
    //! The REAL-BLOCK sibling (Step 3c): a real TRIM block's quadrature runs in REAL arithmetic.
    virtual rsmat_t Matrix(const robs_t* bs, const rvec_t& v) const=0;
    //! \brief \f$\int w_A f\,d^3r\f$ per SITE BLOCK of my quadrature, for a field sampled at MY points -- the
    //! atom-partitioned sibling of \c Integrate.  Default EMPTY: a quadrature with no atomic partition (a
    //! uniform raster) has no basins to integrate over, and the caller must ask rather than assume.
    //! A QUADRATURE question and nothing more (R1.0h): the OBSERVABLE built on it -- the integrated site
    //! moment \f$\mu_A\f$ of \f$\rho_\uparrow-\rho_\downarrow\f$ -- is the spin-native XC term's to
    //! compute and \c ChargeBreakdown's to carry; this engine no longer knows what it is integrating.
    virtual rvec_t SiteIntegrals(const rvec_t& f) const {return rvec_t();}
    //! \brief Pre-warm this engine's per-density caches for \a cd -- the EAGER REFRESH PHASE
    //! (doc/OpenWork.md item **KP**).  \a polarized picks which shape to warm, because the two are
    //! mutually exclusive on one engine (see the cross-invalidation warning on both implementations: an
    //! engine answers \c Rho or \c RhoPol for a run, never both).
    //!
    //! It is expressed as ONE call rather than "the term calls Rho/RhoPol itself" so that the WARMING and
    //! the SHAPE RULE stay in the engine that owns the caches; a term that reached in by name would have
    //! to know the exclusivity rule too.
    virtual void WarmForDensity(const cChargeDensity* cd, bool polarized) const
    {
        if (polarized) { RhoPol(cd, Spin::Up); RhoPol(cd, Spin::Down); }
        else           { Rho(cd); }
    }
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

} //namespace


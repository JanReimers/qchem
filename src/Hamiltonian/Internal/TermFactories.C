// File: Hamiltonian/Internal/TermFactories.C  THE MOLECULAR TERM FACTORIES (V-HAM-COMPOSE step 1, lane S 2026-11-25; design note
// lane-s-V-HAM-COMPOSE-design.md; user 2026-11-25: "What we want to avoid is lots of if statements inside those Hamiltonian constructors.
// Ideally they are a succinct set of easily readable Add(new XXX) calls.").
//
// The concrete wiring the Hamiltonian constructors used to carry inline, moved here unchanged (step 1 is a pure move: no number moves):
//   - MolecularTermContext (was XCMeshShare): the run's XC quadrature and the delta fit on it, ONE per Hamiltonian, shared by the Hartree and
//     XC terms (P2), plus the structure and orbital basis every term factory reads;
//   - MakeCoulombTerm (was MolecularHartreeTerm): THE Hartree route switch, the one place it is decided;
//   - the XC term's nonlinear meshes (MolecularNonlinearMesh / MolecularGradientMesh), the fit-in-the-term's-scope helper, the pin-30
//     admission on built objects, the mesh-client report.
// Factories are allowed to look ugly (CLAUDE.md): the constructors may not.  Step 2 turns the XC side into one MakeXCTerm and the
// constructors into Add lists.
module;
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.Hamiltonian.Internal.TermFactories;
export import qchem.Hamiltonian.Types;
export import qchem.Hamiltonian.Factory;           // st_t, CoulombFitSlot
import qchem.Hamiltonian;                          // rDynamic_HT
import qchem.Mesh.Maker;                           // MeshParams, MeshMaker, SymmetrySource
import qchem.Mesh;                                 // qcMesh::Mesh
import qchem.Mesh.SiteMultipoles;                  // MakeSiteMultipoles
import qchem.ChargeDensity.DensitySampler;         // the samplers the nonlinear meshes are
import qchem.Reporting;                            // BuildTermFit's scope and row
import qchem.Symmetry.Group;                       // RequireForcedWithinG: G->Contains
import qchem.Symmetry.SymOp;                       // InvariantOnly::MeshOps
import qchem.Structure;

export namespace qchem::Hamiltonian
{
//! A FIT BASIS IS BUILT IN ITS TERM'S SCOPE, one level into the caller's context, and files one \c fitbasis row there: a facade's run
//! reads hamiltonian.<term>.fitbasis (user rulings 2026-11-04: only the Hamiltonian knows whether a run needs a fit, so it builds them;
//! the report's nesting is the run sequence).  What the fit reports of its own mesh or grid while it builds lands in the term too, so
//! two fits of one run can never file over each other.  \a term is the report's name for the term the fit will serve.
template <class Fit, class Make> std::shared_ptr<const Fit> BuildTermFit(const std::string& term, const std::string& role, Make make)
{
    qchem::report::Scope here(term);
    std::shared_ptr<const Fit> fit(make());
    qchem::report::Set("fit", fit->Name());   // the resolved fit, on the console (a Scope renders its scalars, not the row)
    qchem::report::Row row("fitbasis");
    qchem::report::Set("role",       role);
    qchem::report::Set("name",       fit->Name());
    qchem::report::Set("nFunctions", (long)fit->GetNumFunctions());
    return fit;
}

//! \brief PIN 30 ADMISSION ON THE BUILT OBJECTS (Group face step 3b; USER RULING 2026-11-21, "I agree with your recommendations": (A)).
//! The built realization and its points each say what they FORCE on the field, as a Group in their own Cartesian frame -- the sampler
//! (\c DensitySamplerBase::Forces: the fit's group, or the periodic A1 space's T) and the XC point set (\c qcMesh::Mesh::Forces, stamped
//! by the structure that built it: a uniform raster T of its cell, the radial-only set O(3) about its site) -- and the structure's own
//! group G must contain both: \f$H_{fit}\subseteq G\f$ and \f$H_{points}\subseteq G\f$.  Every factory asks the cheap pre-flight
//! (\c AdmitVxc, G alone) before it builds, and refuses there; so a cell reaching here was ADMITTED by it, and a refusal here means the
//! pre-flight and the built objects DISAGREE -- a broken invariant naming both answers, a \c logic_error, never a quiet second refusal
//! (controller ruling 2026-11-21, condition (c)).  \a points is null when the fit is its own point set (the plane-wave ball on its
//! raster): its group is the sampler's.
template <class Sampler> const Sampler& RequireForcedWithinG(const Sampler& s, const qcMesh::Mesh* points, const Structure& st)
{
    const Symmetry::GroupPtr G=st.MakeSymmetry(Symmetry::kAdmissionTolerance);
    auto require=[&](const Symmetry::GroupPtr& H, const char* what)
    {
        if (!H) throw std::logic_error(std::string("Vxc builder: the built ")+what+" forces no group (Forces returned null)");
        if (!G->Contains(*H))
            throw std::logic_error(std::string("Vxc builder: the pin-30 pre-flight (Hamiltonian::AdmitVxc, asked of G alone) ADMITTED this "
                                   "cell, but the built ")+what+" forces "+H->Describe()+", which this structure's G ("+G->Describe()+
                                   ") does not contain -- the pre-flight's rows and the built objects disagree (doc/Pins.md pin 30; the drift "
                                   "gate is UTHamiltonian VxcProductSpace)");
    };
    require(s.Forces(), "Vxc realization");
    if (points) require(points->Forces(), "XC point set");
    return s;
}

//! The run's Coulomb fit from its slot (built by \a make on the first ask).
template <class T, class Make> std::shared_ptr<const BasisSet::FIT_CD_ABS<T>> OneCoulombFit(const CoulombFitSlot<T>& slot, Make make)
{
    using fit_t=std::shared_ptr<const BasisSet::FIT_CD_ABS<T>>;
    return slot.GetOrFill([&]{return fit_t(make());});
}

//! THE MOLECULAR / ATOM Vxc DEFAULT: an unresolved \c Auto is the delta basis on the run's own XC mesh (body: Imp/TermFactories.C).
VxcFit MolecularVxcFit(VxcFit fit);
//! The site structure of an atom-centred mesh built from  mp, in the structure's own site order.
std::shared_ptr<const qcMesh::SiteMultipoles> MakeSiteMultipoles(const st_t& st, const qcMesh::Mesh& mesh, const qcMesh::MeshParams& mp, int lMax);
//! THE NONLINEAR MESH'S CLIENTS, listed under its report scope (V-NONLINEAR-MESH).
void ReportMeshClients(bool ppQuadrature, const char* xcClient="Vxc_DFT");
//! The run's nonlinear-integral quadrature at  mp, under the "nonlinearMesh" scope.
BasisSet::FitQuadrature BuildMeshQuadrature(const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry, const qcMesh::MeshParams& mp);

//! \brief \a s's ops WITHOUT the density averaging: a mesh that must be invariant (no deck orientation in the answer) but that no density is
//! averaged on (the PP quadrature).  The factory's wiring, so it lives here.
class InvariantOnly : public qcMesh::SymmetrySource
{
public:
    explicit InvariantOnly(const qcMesh::SymmetrySource& s);
    std::vector<Symmetry::SymOp> MeshOps() const override;
    bool AveragesDensityOverOrbits() const override;
private:
    const qcMesh::SymmetrySource& itsS;
};


//! \brief THE XC QUADRATURE AND THE δ FIT BASIS ON IT, ONE PER HAMILTONIAN (lane H 2026-11-15, the basis-on-grid survey P2,
//! controller ruling).  A molecular run whose Hartree term is a δ-fit route (\c CoulombFit::Delta / \c Multipole) on the XC mesh's own
//! parameters used to build a SECOND quadrature and a second δ fit beside the Vxc term's, identical, so every orbital block's
//! \f$\Phi_{gi}=\chi_i(r_g)\f$ table (npts x n, the largest grid holder) was held twice.  The factory now builds them once and hands
//! the same objects to both terms: a pairing decided here, where the concrete wiring lives; neither term caches anything new.  Each is
//! built on the first ask, inside the asking term's report scope, so the run report shows which term built it.
class MolecularTermContext
{
public:
    MolecularTermContext(const st_t& st, const rbs_t* bs, const qcMesh::MeshParams& mp, const qcMesh::MeshMaker& meshMaker, const qcMesh::SymmetrySource& symmetry);
    const qcMesh::MeshParams& Params() const {return itsMp;}
    const qcMesh::MeshMaker&      MeshMaker() const {return itsMeshMaker;}
    const qcMesh::SymmetrySource& Symmetry () const {return itsSymmetry;}
    const BasisSet::FitQuadrature& Quadrature();
    //! The δ fit basis on \c Quadrature(); \a built tells the caller whether this ask built it (else another term did).
    std::shared_ptr<const BasisSet::rFIT_SF_ABS> Delta(bool& built);

    //! \name THE PER-ORBIT XC ROUTE (OpenWork §3, the cost lever V-FOLD-NEUTRAL did not buy; lane H 2026-11-24)
    //! On a FOLDED quadrature (an averaging run: the SALC-imposed molecule) the δ XC route samples ONE POINT PER ORBIT at the orbit's summed
    //! weight (\c FoldedMesh::MakeRepresentativeMesh).  Exact, not an approximation: every integrand the δ route forms within a block -- ρ,
    //! \f$v\,\chi_i\chi_j\f$, \f$|\nabla\rho|^2\f$, \f$\nabla\rho\cdot\nabla(\chi_i\chi_j)\f$ -- is invariant under the fold's group
    //! when each block is a 1-D real irrep of it, and the δ fit function at a representative with weight \f$W_o\f$ IS the orbit sum
    //! \f$\sum_R\delta(R\,r_o)\f$.  So no star-average is left to do (the quadrature handed over carries no fold), and a GGA needs no
    //! rotation of \f$\nabla\rho\f$.  ONLY the δ route: a Gaussian or multipole Vxc fit projects onto single fit functions, which are not
    //! invariant, so those keep the full folded quadrature (\c Quadrature()).  The δ-fit Hartree samples on the SAME representatives
    //! (\c XCDelta(), shared): its operator couples every point, but an invariant charge has an invariant potential, so the operator scatters
    //! the orbits' charges to the full mesh, solves there and reads v back at the representatives (\c qcMesh::MakeFoldedCoulomb).
    //!@{
    //! The quadrature the δ XC route samples on: the representatives when the run folds, else \c Quadrature().
    const BasisSet::FitQuadrature& XCQuadrature();
    //! The δ fit basis on \c XCQuadrature(): the shared \c Delta() on an unfolded run, its own on the representatives on a folded one.
    std::shared_ptr<const BasisSet::rFIT_SF_ABS> XCDelta(bool& built);
    //!@}
    //! \brief THROWS unless every orbital block is a 1-D irrep of the fold's group: the per-orbit route's precondition.  The policy layer
    //! folds only a SALC-imposed run, so a failure here is a broken invariant.  The block's own \c Symmetry::Representation answers it:
    //! its image of a generic symmetric matrix under every op must be the matrix itself (a SALC block is the identity on its stabiliser,
    //! \f$\chi(g)^2=1\f$; an AO block permutes and rotates; a block that does not know the op throws).
    void RequireInvariantPairProducts() const;
private:

    st_t                                         itsSt;
    const rbs_t*                                 itsBs;
    qcMesh::MeshParams                           itsMp;
    const qcMesh::MeshMaker&                     itsMeshMaker;   //!< the factory's caller's, alive for the Hamiltonian's construction
    const qcMesh::SymmetrySource&                itsSymmetry;    //!< likewise
    BasisSet::FitQuadrature                      itsQuad;
    std::shared_ptr<const BasisSet::rFIT_SF_ABS> itsDelta;
    BasisSet::FitQuadrature                      itsRepQuad;    //!< the representatives (folded runs only; empty otherwise)
    std::shared_ptr<const BasisSet::rFIT_SF_ABS> itsRepDelta;   //!< the δ fit on them
};


//! THE MOLECULAR / ATOM NONLINEAR INTEGRATION MESH: the sampler the Hamiltonian holds and the ONE XC term runs on (V-NONLINEAR-MESH).
using rmesh_t=std::shared_ptr<const ChargeDensity::tDensitySampler<double>>;
rmesh_t MolecularNonlinearMesh(const st_t& st, const rbs_t* bs, MolecularTermContext& share, XCRoute xc, int lMax);
//! THE MOLECULAR DENSITY-GRADIENT MESH (GGA stage 2), on the delta Vxc fit only.
using rgmesh_t=std::shared_ptr<const ChargeDensity::tDensityGradientSampler<double>>;
rgmesh_t MolecularGradientMesh(const st_t& st, MolecularTermContext& share, XCRoute xc);
//! THE MOLECULAR / ATOM HARTREE TERM: the ONE place the Hartree route is switched (Auxiliary Gaussian fit | delta: radial Poisson, direct
//! sum, analytic (ij|delta) | multipole).   coulombFit: the run's Coulomb-fit slot, filled on the first build, read on every later one.
rDynamic_HT* MakeCoulombTerm(const st_t& st, const rbs_t* bs, MolecularTermContext& share, const CoulombRoute& coulomb,
                             const CoulombFitSlot<double>& coulombFit);
}

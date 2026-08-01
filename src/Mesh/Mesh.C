// File: Mesh.C  Clean-room replacement for qcMesh.  The geometry-free Mesh value type.
//
// A Mesh is a concrete VALUE type = quadrature points + weights, stored SEPARATELY (SoA).
// Points are streamed by evaluators (phi(r)); weights by integrators.  They live in totally
// different algorithms, so they get totally separate arrays.  No polymorphism, no Clone.
//
// Everything lives in namespace qchem::qcMesh so the new library coexists with the old qchem.Mesh
// (which exports a global `class Mesh`) during the migration; the namespace goes away once the
// old library is deleted.
module;
#include <utility>
#include <cassert>
#include <string>
export module qchem.Mesh;
export import qchem.Types;

export namespace qchem::qcMesh
{

//! \brief A quadrature mesh: points r_i and weights w_i stored as separate arrays (SoA).
//!
//! \f$\intop f\,d^3r\approx\sum_i w_i f(r_i)\f$.  This is a plain value type with no virtuals:
//! builders \c Append points, consumers read \c Points() (in op(r) evaluators) and \c Weights()
//! (in the integration algorithms) independently.
class Mesh
{
public:
    Mesh() = default;
    //! Build from ready-made SoA arrays (the efficient path when the size is known up front,
    //! e.g. ProductMesh = nRadial*nAngular).  Sizes must match.
    Mesh(rvec3vec_t r, rvec_t w) : itsR(std::move(r)), itsW(std::move(w)) {}

    const rvec3vec_t& Points () const {return itsR;} //!< r_i, for the phi(r) evaluators.
    const rvec_t&     Weights() const {return itsW;} //!< w_i, for the integrators.
    size_t            size   () const {return itsW.size();}

    //! Append one (point,weight) pair.  Builders (the incremental Becke mesh) push here.
    void Append(const rvec3_t& r, double w);
    //! Translate every point: r_i += o.  Used to place a single-center mesh at an atom.
    void ShiftOrigin(const rvec3_t& o);

private:
    rvec3vec_t itsR;
    rvec_t     itsW;
};

//! \brief Radial mesh family.  See the per-class transplanted formulae in Internal/.
enum class RadialKind  {MHL, Log, Linear};
//! \brief Angular mesh family.  All schemes are normalised so \f$\sum_i w_i = 4\pi\f$.
enum class AngularKind {Gauss, GaussLegendre, EulerMaclaren};
//! \brief Which quadrature a periodic UnitCell builds: the uniform FFT-compatible midpoint grid
//! (Hartree/PP -- resolution from \c eCut / \c nUniform), or the atom-centred periodic Becke
//! fuzzy-Voronoi grid (dense radial near each nucleus, cheap diffuse tails -- the XC quadrature;
//! uses the \c radial / \c nAngular / \c beckeOrder knobs).  See doc/GPWPlan1.md "Becke XC grid".
//! \brief Which quadrature a UnitCell builds for a lattice mesh.  \c Auto = "the caller did not choose":
//! a POLICY layer that knows the run context resolves it (the GPW driver picks Becke unless the run is
//! BZ-reduced -- doc/GPWPlan1.md "Becke as the DEFAULT"); any consumer reached WITHOUT resolution treats
//! \c Auto as the historical \c Uniform (every test below compares \c ==Becke, so Auto falls through safely).
enum class UnitCellKind {Uniform, Becke, Auto};

//! \brief Becke's iterated smoothing polynomial mapped to the cell cutoff:
//! \f$ s(\mu) = \tfrac12(1-f_k(\mu)),\; f(\mu)=\tfrac12(3\mu-\mu^3) \f$ applied \f$k+1\f$ times
//! (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)).  \f$ s:[-1,1]\to[1,0] \f$, with \f$s(\mp1)\f$
//! pinned at 1/0 and the transition sharpened by each iteration.  Shared by the molecular
//! (MakeMolecularMesh) and periodic (UnitCell) fuzzy-Voronoi partitions.
double BeckeCutoff(double mu, int k);

//! \brief Typed, fully-defaulted mesh parameters.  Set only the knobs you actually use
//! (C++20 designated initializers); no field is required-but-unused.
struct MeshParams
{
    RadialKind  radial    = RadialKind::MHL;    int    nRadial   = 30;
    int         mhl_m     = 2;                   double mhl_alpha = 1.0;     //!< MHL only.
    double      logStart  = 1.0e-4;             double logStop   = 50.0;     //!< Log only.
    AngularKind angular   = AngularKind::Gauss;  int    nAngular  = 12;      //!< Gauss: #dirs; GL/EM: L.
    int         em_m      = 2;                                               //!< EulerMaclaren only (1..3).
    int         beckeOrder= 3;   //!< Becke fuzzy-Voronoi smoothing iterations (molecular mesh only).
    int         nUniform  = 20;  //!< Uniform periodic real-space grid: points per cell axis (lattice mesh only; \f$n^3\f$ total). Manual fallback when \c eCut<=0.
    double      eCut      = 0.0;  //!< Real-space integration-mesh energy cutoff (a.u.). If >0, the uniform lattice mesh DERIVES its \c nUniform from the Nyquist bound \f$n\gtrsim 2a\sqrt{2E_{cut}}/\pi\f$ (\f$a\f$=longest cell edge; the \f$\times2\f$ is the density bandwidth), and \c nUniform is ignored. 0=use the manual \c nUniform.
    double      relCutoff = 1.0;  //!< Fit-grid density multiplier (CP2K \c REL_CUTOFF): the fit basis scales its \f$E_{cut}\f$ by this. 1=wavefunction bandwidth (LDA); GGA wants >1. Set by the Hamiltonian from the functional's \c GridCutoffFactor().
    UnitCellKind cellKind = UnitCellKind::Uniform;  //!< Lattice mesh only: which quadrature a UnitCell builds (uniform midpoint grid vs periodic Becke).

    //! \brief Compact, deterministic identity string for these parameters.  Two MeshParams give the
    //! same ID() iff they build the same quadrature, so it is the cache key for any mesh-quadrature
    //! quantity that must NOT be shared across different meshes (e.g. a fit basis's numerical Norm --
    //! same basis, different mesh => different normalisation).  All fields are folded in; over-keying
    //! on an inactive field (e.g. logStart when radial=MHL) is harmless.
    std::string ID() const
    {
        using std::to_string;
        return "M{r"  + to_string(static_cast<int>(radial))  + ",nr" + to_string(nRadial)
             + ",m"   + to_string(mhl_m)    + ",a"  + to_string(mhl_alpha)
             + ",ls"  + to_string(logStart) + ",le" + to_string(logStop)
             + ",ang" + to_string(static_cast<int>(angular)) + ",na" + to_string(nAngular)
             + ",em"  + to_string(em_m)     + ",bo" + to_string(beckeOrder)
             + ",nu"  + to_string(nUniform) + ",ec" + to_string(eCut)
             + ",rc" + to_string(relCutoff) + ",ck" + to_string(static_cast<int>(cellKind)) + "}";
    }
};

} //export namespace qchem::qcMesh

//-----------------------------------------------------------------------------------------------
namespace qchem::qcMesh
{

// Append grows the SoA arrays by one.  blaze resize copies, so this is a setup-path convenience
// (the incremental Becke builder), NOT for hot loops -- when the size is known, use the
// from-arrays constructor instead.
void Mesh::Append(const rvec3_t& r, double w)
{
    size_t n=itsW.size();
    assert(itsR.size()==n);
    itsR.resize(n+1,true);
    itsW.resize(n+1,true);
    itsR[n]=r;
    itsW[n]=w;
}

void Mesh::ShiftOrigin(const rvec3_t& o)
{
    for (size_t i=0; i<itsR.size(); i++) itsR[i]+=o;
}

double BeckeCutoff(double mu, int k)
{
    for (int i=k; i>=0; i--) mu=0.5*(3*mu - mu*mu*mu);
    return 0.5*(1.0-mu);
}

} //namespace qchem::qcMesh

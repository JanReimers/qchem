// File: Mesh1.C  Clean-room replacement for qcMesh.  The geometry-free Mesh value type.
//
// A Mesh is a concrete VALUE type = quadrature points + weights, stored SEPARATELY (SoA).
// Points are streamed by evaluators (phi(r)); weights by integrators.  They live in totally
// different algorithms, so they get totally separate arrays.  No polymorphism, no Clone.
//
// Everything lives in namespace qcMesh1 so the new library coexists with the old qchem.Mesh
// (which exports a global `class Mesh`) during the migration; the namespace goes away once the
// old library is deleted.
module;
#include <utility>
#include <cassert>
export module qchem.Mesh1;
export import qchem.Types;

export namespace qcMesh1
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
};

} //export namespace qcMesh1

//-----------------------------------------------------------------------------------------------
namespace qcMesh1
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

} //namespace qcMesh1

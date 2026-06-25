// File: Fields.C  Pointwise real-space field interfaces.
//
// These are deliberately INDEPENDENT of Mesh.  In the old qcMesh, ScalarFunction/VectorFunction
// carried a virtual operator()(const Mesh&) -- an ISP violation that dragged Mesh into every
// pointwise field and forced a world rebuild whenever Mesh.C changed.  Here a field knows only
// how to answer at a single point r; the free-function quadrature (qchem.Mesh1.Quadrature) is
// what streams a Mesh through it.
module;
export module qchem.Mesh1.Fields;
export import qchem.Types;

export namespace qcMesh1
{

//! \brief A scalar field: rho(r), vxc(r), -Z/|r-R|, 1/r, ...  Geometry/physics live in the caller.
template <class T> class ScalarField
{
public:
    virtual ~ScalarField() = default;
    virtual T         operator()(const rvec3_t&) const = 0;   //!< f(r)
    virtual vec3_t<T> Gradient  (const rvec3_t&) const = 0;   //!< grad f(r)
};

//! \brief A basis as a vector field: [phi_i(r)].  size() is the number of basis functions
//! (needed by the quadrature to size its matrices -- this is the one addition over the design sketch).
template <class T> class BasisField
{
public:
    virtual ~BasisField() = default;
    virtual size_t       size()                     const = 0;   //!< number of basis functions
    virtual vec_t<T>     operator()(const rvec3_t&)  const = 0;   //!< [ phi_i(r) ]
    virtual vec3vec_t<T> Gradient  (const rvec3_t&)  const = 0;   //!< [ grad phi_i(r) ]
};

} //export namespace qcMesh1

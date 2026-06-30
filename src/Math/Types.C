// File: Types.C  Some common types used everywhere in the code
module;
#include <blaze/Math.h>
#include <complex>
#include <ranges>
#include <type_traits>
export module qchem.Types;
export import qchem.Vector3D;
export import qchem.Matrix3D;   // fixed 3x3 (brings Vector3D too)
export import qchem.Matrix2D;   // brings in Vector2D too (fixed 2x2 / 2-vector)

export
{
// We need size_t everywhere.  Clang takes all three options.
// export using size_t = decltype(sizeof 0); //gcc-15-1 accepts this.
using std::size_t; //gcc-15-1 rejects this.
// export typedef std::size_t size_t; //gcc-15-1 rejects this.

using dcmplx=std::complex<double>;

template <typename T> using vec3_t = qchem::Vector3D<T>;
template <typename T> using mat3d_t = qchem::Matrix3D<T>;
template <typename T> using vec2_t = qchem::Vector2D<T>;
template <typename T> using mat2d_t = qchem::Matrix2D<T>;
template <typename T> using  mat_t = blaze::DynamicMatrix<T,blaze::columnMajor>;
template <typename T> using smat_t = blaze::SymmetricMatrix<mat_t<T>>;
template <typename T> using hmat_t = std::conditional_t<
    std::is_floating_point_v<T>,
    blaze::SymmetricMatrix<mat_t<T>>,
    blaze::HermitianMatrix<mat_t<T>>>;
template <typename T> using  vec_t = blaze::DynamicVector<T>;
template <typename T> using  row_t = blaze::DynamicVector<T,blaze::rowVector>;
template <typename T> using  col_t = blaze::DynamicVector<T,blaze::columnVector>;
template <typename T> using vec3vec_t = vec_t<vec3_t<T>>;
template <typename T> using vec3mat_t = mat_t<vec3_t<T>>;

using  rvec_t= vec_t<double>;
using  ivec_t= vec_t<int>;
using  rrow_t= row_t<double>;
using  rcol_t= col_t<double>;
using  rmat_t= mat_t<double>;
using rsmat_t=smat_t<double>;
using rhmat_t=hmat_t<double>;
using rvec3_t=vec3_t<double>;
using rmat3d_t=mat3d_t<double>;    // fixed 3x3 matrix (e.g. inertia tensor, unit cell, rotations)
using rvec3vec_t=vec3vec_t<double>;
using ivec3_t = vec3_t<int>;
using rvec2d_t=vec2_t<double>;     // fixed 2-vector (e.g. LAPW matching coefficients)
using rmat2d_t=mat2d_t<double>;    // fixed 2x2 matrix (e.g. LAPW {u,udot}-basis blocks)

// Complex (dcmplx) counterparts of the real aliases above, for plane-wave / Bloch work.
// Hermitian (chmat_t), not symmetric: a complex-symmetric matrix is rarely physical, whereas
// overlap / kinetic / potential blocks are Hermitian.  (chmat_t mirrors rhmat_t.)
using chmat_t=hmat_t<dcmplx>;
using cvec_t=vec_t<dcmplx>;
using cvec3vec_t=vec3vec_t<dcmplx>;

// Used for angular intgrals Ak arrays.
typedef blaze::StaticVector<double,11> rvec11_t; //la+lc+1<=11 support up to i orbitals.  Good luck finding a stable nucleus!
typedef std::ranges::iota_view<size_t,size_t> iv_t; //For range based loops
 

}
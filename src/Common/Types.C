// File: Types.C  Some common types used everywhere in the code
module;
#include "blaze/Math.h" 
#include <complex>
#include <ranges>
export module qchem.Types;
export import qchem.Vector3D;

export
{
// We need size_t everywhere.  Clang takes all three options.
// export using size_t = decltype(sizeof 0); //gcc-15-1 accepts this.
using std::size_t; //gcc-15-1 rejects this.
// export typedef std::size_t size_t; //gcc-15-1 rejects this.

using dcmplx=std::complex<double>;

template <typename T> using vec3_t = Vector3D<T>;
template <typename T> using  mat_t = blaze::DynamicMatrix<T,blaze::columnMajor>;
template <typename T> using smat_t = blaze::SymmetricMatrix<mat_t<T>>;
template <typename T> using  vec_t = blaze::DynamicVector<T>;
template <typename T> using  row_t = blaze::DynamicVector<T,blaze::rowVector>;
template <typename T> using  col_t = blaze::DynamicVector<T,blaze::columnVector>;
template <typename T> using vec3vec_t = vec_t<vec3_t<T>>;
template <typename T> using vec3mat_t = mat_t<vec3_t<T>>;

using  rvec_t= vec_t<double>;
using  rrow_t= row_t<double>;
using  rcol_t= col_t<double>;
using  rmat_t= mat_t<double>;
using rsmat_t=smat_t<double>;
using rvec3_t=vec3_t<double>;
using rvec3vec_t=vec3vec_t<double>;
using ivec3_t = vec3_t<int>;

typedef std::ranges::iota_view<size_t,size_t> iv_t; //For range based loops


}
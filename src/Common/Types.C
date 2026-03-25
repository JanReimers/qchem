// File: Types.C  Some common types used everywhere in the code
module;
#include "blaze/Math.h" 
#include <complex>
#include <ranges>
export module qchem.Types;
export import oml.Vector3D;
export import oml.Vector;
export import oml.Matrix;
export import oml.SMatrix;

export
{
// We need size_t everywhere.  Clang takes all three options.
// export using size_t = decltype(sizeof 0); //gcc-15-1 accepts this.
using std::size_t; //gcc-15-1 rejects this.
// export typedef std::size_t size_t; //gcc-15-1 rejects this.

using dcmplx=std::complex<double>;

typedef Vector3D<double> RVec3;
typedef Vector3D<dcmplx> CVec3;
typedef Vector3D<int>    IVec3;
typedef Vector3D<size_t> UVec3;

typedef Vector<double>   RVec;

template <typename T> using  mat_t = blaze::DynamicMatrix<T,blaze::columnMajor>;
template <typename T> using smat_t = blaze::SymmetricMatrix<mat_t<T>>;
template <typename T> using  vec_t = blaze::DynamicVector<T>;
template <typename T> using  row_t = blaze::DynamicVector<T,blaze::rowVector>;
template <typename T> using  col_t = blaze::DynamicVector<T,blaze::columnVector>;

using  rvec_t= vec_t<double>;
using  rrow_t= row_t<double>;
using  rcol_t= col_t<double>;
using  rmat_t= mat_t<double>;
using rsmat_t=smat_t<double>;

typedef std::ranges::iota_view<size_t,size_t> iv_t; //For range based loops


}
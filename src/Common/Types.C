// File: Types.C  Some common types used everywhere in the code
module;
#include "blaze/Math.h" 
#include <complex>
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

using  rvec_t= vec_t<double>;
using  rmat_t= mat_t<double>;
using rsmat_t=smat_t<double>;

template <typename T> smat_t<T> convert(const SMatrix<T>& S)
    {
        size_t N=S.GetNumRows();
        smat_t<T> bS(N);
        for (auto i:S.rows())
                for (auto j:S.cols(i))
                    bS(i-1,j-1)=S(i,j);
        return bS;
    }
    template <typename T> mat_t<T> convert(const Matrix<T>& M)
    {
        size_t N=M.GetNumRows();
        mat_t<T> bM(N,N);
        for (auto i:M.rows())
                for (auto j:M.cols())
                    bM(i-1,j-1)=M(i,j);
        return bM;
    }
    template <typename T> SMatrix<T> convert(const smat_t<T>& bS)
    {
        SMatrix<T> S(bS.rows());
        for (auto i:S.rows())
                for (auto j:S.cols(i))
                    S(i,j)=bS(i-1,j-1);
        return S;
    }
    template <typename T> Matrix<T> convert(const mat_t<T>& bM)
    {
        Matrix<T> M(bM.rows(),bM.columns());
        for (auto i:M.rows())
                for (auto j:M.cols())
                    M(i,j)=bM(i-1,j-1);
        return M;
    }
    template <typename T>  vec_t<T> convert(const Vector<T>& V)
    {
        vec_t<T> bV(V.size());
        for (auto i:bV.indices())
            V[i-1]=bV(i);
        return bV;
    }
    template <typename T> Vector<T> convert(const vec_t<T>& bV)
    {
        Vector<double> V(bV.size());
        for (auto i:V.indices())
                    V(i)=bV[i-1];
        return V;
    }

}
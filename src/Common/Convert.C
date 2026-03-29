// File: Common/Convert.C
module;
#include <vector>
#include <valarray>
#include "blaze/Math.h" 
export module qchem.Conversions;
export import qchem.Types;
import oml.Vector;
import oml.Matrix;
import oml.SMatrix;

export 
{

template <class T> std::valarray<T> to_valarray(const std::vector<T>& v)
{
    std::valarray<T> ret(v.size());
    size_t n=0;
    for (auto iv:v) ret[n++]=iv;
    return ret;
}

template <class T> Vector<T> to_omlVector(const std::valarray<T>& v)
{
    Vector<T> ret(v.size());
    size_t n=0;
    for (auto i:v) ret(++n)=i;
    return ret;
}

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
        for (auto i:V.indices())
            bV[i-1]=V(i);
        return bV;
    }
    template <typename T>  vec_t<T> convert1(const std::valarray<T>& V)
    {
        vec_t<T> bV(V.size());
        for (auto i:iv_t(0,V.size()))
            bV[i]=V[i];
        return bV;
    }
    template <typename T> Vector<T> convert(const vec_t<T>& bV)
    {
        Vector<double> V(bV.size());
        for (auto i:V.indices())
                    V(i)=bV[i-1];
        return V;
    }

    template <typename T> smat_t<T> zero(size_t N)
    {
        smat_t<T> z(N);
        for ( size_t j=0; j<z.columns(); ++j) column(z,j)=0.0;
        return z;
    }

    template <typename T> double fnorm(const smat_t<T>& S)
    {
        return sqrt(sum(S%conj(S)));
    }


}


// File: Common/Convert.C
module;
#include <vector>
#include <valarray>
#include "blaze/Math.h" 
export module qchem.Conversions;
export import qchem.Types;

export 
{

    template <class T> std::valarray<T> to_valarray(const std::vector<T>& v)
    {
        std::valarray<T> ret(v.size());
        size_t n=0;
        for (auto iv:v) ret[n++]=iv;
        return ret;
    }

    template <typename T>  vec_t<T> convert1(const std::valarray<T>& V)
    {
        vec_t<T> bV(V.size());
        for (auto i:iv_t(0,V.size()))
            bV[i]=V[i];
        return bV;
    }

    template <typename T> smat_t<T> zero(size_t N)
    {
        smat_t<T> z(N);
        for ( size_t j=0; j<z.columns(); ++j) column(z,j)=0.0;
        return z;
    }

    

}


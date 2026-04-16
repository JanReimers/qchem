// File: Common/Convert.C
module;
#include <blaze/math/SymmetricMatrix.h>
export module qchem.Blaze;
export import qchem.Types;

export 
{
    template <typename T> smat_t<T> zero(size_t N)
    {
        smat_t<T> z(N);
        for ( size_t j=0; j<z.columns(); ++j) column(z,j)=0.0;
        return z;
    }

}


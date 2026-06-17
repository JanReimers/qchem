// File: Common/Convert.C
module;
#include "blaze/math/dense/DenseVector.h"
#include <blaze/Math.h>
export module qchem.Blaze;
export import qchem.Types;

export 
{
    using blaze::DenseIterator;
    using blaze::operator==;

    using blaze::operator+;
    using blaze::operator/;
    using blaze::operator*;
    using blaze::operator-;

    using blaze::operator+=;
    using blaze::operator-=;
    using blaze::operator*=;
    using blaze::operator/=;

    using blaze::operator%;

    using blaze::operator<<;

}


export namespace blazem
{
    template <typename T> smat_t<T> zero(size_t N)
    {
        smat_t<T> z(N);
        for ( size_t j=0; j<z.columns(); ++j) column(z,j)=0.0;
        return z;
    }

    
    using blaze::DynamicVector;
    using blaze::Subvector;
    using blaze::subvector;
    using blaze::submatrix;
    using blaze::exp;
    using blaze::acos;
    using blaze::abs;
    using blaze::max;
    using blaze::min;
    using blaze::size;
    using blaze::isnan;
    using blaze::sqrt;
    using blaze::column;
    using blaze::sum;
    using blaze::norm;
    using blaze::linspace;
}


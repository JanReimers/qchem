// File: Common/Convert.C
module;
#include <blaze/Math.h>
export module qchem.Blaze;
export import qchem.Types;
//
//  Export some blaze functions into a new namespace blazem in order get them into a module BMI file.
//  This is basically a way of making a pre-compiled header for some blaze function we need.
//
export namespace blazem
{
    template <typename T> smat_t<T> zero(size_t N)
    {
        smat_t<T> z(N);
        for ( size_t j=0; j<z.columns(); ++j) column(z,j)=0.0;
        return z;
    }
   
    using blaze::DynamicVector;
    using blaze::DiagonalMatrix;
    using blaze::UpperMatrix;
    using blaze::LowerMatrix;
    using blaze::Subvector;
    using blaze::eigen;
    using blaze::svd;
    using blaze::inv;
    using blaze::solve;
    using blaze::clear;
    using blaze::subvector;
    using blaze::submatrix;
    using blaze::trans;
    using blaze::exp;
    using blaze::acos;
    using blaze::abs;
    using blaze::max;
    using blaze::min;
    using blaze::size;
    using blaze::isnan;
    using blaze::sqrt;
    using blaze::column;
    using blaze::row;
    using blaze::sum;
    using blaze::norm;
    using blaze::outer;
    using blaze::dot;
    using blaze::linspace;
    using blaze::diagonal;
    using blaze::isZero;
    using blaze::isSquare;
    using blaze::conj;
    using blaze::potrf;
    using blaze::trtri;
    
}

// Export a number of overloaded operators.  These have to go into the global namespace.
export 
{
    using blaze::DenseIterator;
    using blaze::operator==;
    using blaze::operator!=;
    // using blaze::operator<;
    // using blaze::operator<=;
    // using blaze::operator>;
    // using blaze::operator>=;

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

    using blaze::begin;
    using blaze::end;
}


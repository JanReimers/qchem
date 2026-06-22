// File: Common/Convert.C
module;
#include <blaze/Math.h>
#include <utility>
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

    // Hermitian counterpart of zero(): an NxN zeroed hmat_t (SymmetricMatrix for real T,
    // HermitianMatrix for complex T).  A zero matrix is trivially Hermitian, so the per-column
    // scalar zeroing keeps the adaptor's invariant.
    template <typename T> hmat_t<T> zeroH(size_t N)
    {
        hmat_t<T> z(N);
        for ( size_t j=0; j<z.columns(); ++j) column(z,j)=0.0;
        return z;
    }

    //
    //  Accumulate a blaze vector when the final length is not known up front.  push_back() grows by
    //  doubling, so N appends cost O(N) total -- a naive resize-by-1 each time would be O(N^2) since
    //  blaze resize copies the whole vector.  take() shrinks to the logical length and hands over the
    //  vector (the builder is empty afterwards).  If you DO know the length, just size the vector once
    //  (`vec_t<T> v(N); v[i]=...;`) -- no builder needed.
    //  NOTE: this is a convenience for setup/parse paths.  It is NOT for hot code.
    //
    template <typename T> class VecBuilder
    {
    public:
        explicit VecBuilder(size_t reserve=8) : itsVec(reserve<1?1:reserve), itsN(0) {}
        void     Append(const T& x)
        {
            if (itsN>=itsVec.size()) itsVec.resize(itsVec.size()*2,true);
            itsVec[itsN++]=x;
        }
        vec_t<T> take()       { itsVec.resize(itsN,true); return std::move(itsVec); }
        size_t   size() const { return itsN; }
    private:
        vec_t<T> itsVec;
        size_t   itsN;
    };
    template <typename T> inline void push_back(VecBuilder<T>& b, const T& x) { b.Append(x); }


    using blaze::DynamicVector;
    using blaze::DiagonalMatrix;
    using blaze::UpperMatrix;
    using blaze::LowerMatrix;
    using blaze::HermitianMatrix;
    using blaze::ctrans;
    using blaze::IdentityMatrix;
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
    using blaze::real;
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


// File: BlazeCore.C  The ONE owner of Blaze's headers (module qchem.Blaze.Core).
//
// Blaze is textual (header-only expression templates), so whichever module interface #includes it absorbs it into its BMI.
// Before this module two interfaces did (qchem.Types and qchem.Blaze), so every importer of both saw two copies of each Blaze
// partial specialization, and clang 23's reduced BMIs could not merge them (MultTraitEval2, SolveTraitEval2 ... ambiguous).
// Here the headers sit in exactly one BMI, and what the rest of the code uses is EXPORTED by name -- the same pattern as
// libc++'s std.cppm (`export namespace std { using std::vector; }`) -- so a client never depends on reaching a declaration
// through another module's global-module fragment.
module;
#include <blaze/Math.h>
export module qchem.Blaze.Core;

export namespace blaze
{
    // Types and tags (qchem.Types' aliases, qchem.Blaze's helpers).
    using blaze::DynamicMatrix;
    using blaze::DynamicVector;
    using blaze::StaticVector;
    using blaze::SymmetricMatrix;
    using blaze::HermitianMatrix;
    using blaze::DiagonalMatrix;
    using blaze::UpperMatrix;
    using blaze::LowerMatrix;
    using blaze::IdentityMatrix;
    using blaze::Subvector;
    using blaze::DenseIterator;
    using blaze::blas_int_t;
    using blaze::columnMajor;    // storage-order / transpose flags: `inline constexpr` in the fork (topic (d)), so exportable
    using blaze::rowMajor;
    using blaze::columnVector;
    using blaze::rowVector;
    // Functions (qchem.Blaze re-exports them into qchem::blazem).
    using blaze::ctrans;
    using blaze::eigen;
    using blaze::heevx;
    using blaze::syevx;
    using blaze::declherm;
    using blaze::declsym;
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
    using blaze::isSymmetric;
    using blaze::isHermitian;
    // Blaze's own free functions, the ones its templates call UNQUALIFIED (ADL) at instantiation.  A full BMI keeps every
    // global-module-fragment declaration anyway; a REDUCED BMI keeps only what the module names, so this list is what lets
    // reduced BMIs build (measured: without it 51, then 39 TUs fail).  Add to it when Blaze grows a predicate qchem reaches.
    using blaze::isIntact;
    using blaze::isDefault;
    using blaze::isSame;
    using blaze::isDivisor;
    using blaze::isOne;
    using blaze::isStrictlyLower;
    using blaze::isStrictlyUpper;
    using blaze::isDiagonal;
    using blaze::isReal;
    using blaze::isUpper;
    using blaze::isLower;
    using blaze::isUniform;
    using blaze::isUniUpper;
    using blaze::isUniLower;
    using blaze::isIdentity;
    using blaze::isPositiveDefinite;
    using blaze::isEmpty;
    using blaze::isinf;
    using blaze::isfinite;
    using blaze::trace;
    using blaze::det;
    using blaze::rank;
    using blaze::lu;
    using blaze::llh;
    using blaze::qr;
    using blaze::rq;
    using blaze::ql;
    using blaze::lq;
    using blaze::reset;
    using blaze::resize;
    using blaze::shrinkToFit;
    using blaze::swap;
    using blaze::transpose;
    using blaze::ctranspose;
    using blaze::invert;
    using blaze::serial;
    using blaze::capacity;
    using blaze::rows;
    using blaze::columns;
    using blaze::conj;
    using blaze::real;
    using blaze::potrf;
    using blaze::pstrf;
    using blaze::trtri;
    using blaze::begin;
    using blaze::end;
    // Operators (qchem.Blaze re-exports them into the global namespace).
    using blaze::operator==;
    using blaze::operator!=;
    using blaze::operator+;
    using blaze::operator-;
    using blaze::operator*;
    using blaze::operator/;
    using blaze::operator+=;
    using blaze::operator-=;
    using blaze::operator*=;
    using blaze::operator/=;
    using blaze::operator%;
    using blaze::operator<<;
}

// File: Atom/ml/BSpline_BF.C  B-Spline basis function.
module;

module qchem.BasisSet.Atom.Internal.ml.BSplineBS;

namespace Atom_ml
{
namespace BSpline
{

template <size_t K> BasisFunction<K>::BasisFunction(const spline_t& sp, int l, int _ml, double norm)
    : Base(sp,l,norm)
    , ml(_ml)
{
    
};


template <size_t K> ::Real_BF* BasisFunction<K>::Clone() const
{
    return new  BasisFunction<K>(*this);
}

#define INSTANCEk(k) template class BasisFunction<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace

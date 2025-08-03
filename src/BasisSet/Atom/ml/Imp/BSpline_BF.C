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



#define INSTANCEk(k) template class BasisFunction<k>;
#include "../../radial/BSpline/Instance.hpp"

}} //namespace

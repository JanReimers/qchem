// File: Atom/ml/Gaussian_BF.C r^l exp(-a*r^2) type Gaussian basis function.

#include "Imp/BasisSet/Atom/ml/Gaussian_BF.H"
#include <cassert>

namespace Atom_ml
{
namespace Gaussian
{

BasisFunction::BasisFunction(double e,int n, int l, int _ml, double norm)
: Base(e,l,norm)
, ml(_ml)
{};

bool BasisFunction::operator==(const ::BasisFunction& bf) const
{
    const BasisFunction& sgbf = dynamic_cast<const BasisFunction&>(bf);
    assert(&sgbf);
    return  Base::operator==(bf) && (ml==sgbf.ml);
}

BasisFunction* BasisFunction::Clone() const
{
    return new  BasisFunction(*this);
}

}}//namespace
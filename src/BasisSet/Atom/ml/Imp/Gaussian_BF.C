// File: Atom/ml/Gaussian_BF.C r^l exp(-a*r^2) type Gaussian basis function.
module;
#include <cassert>
module qchem.BasisSet.Atom.Internal.ml.GaussianBS;

namespace Atom_ml
{
namespace Gaussian
{

BasisFunction::BasisFunction(double e,int n, int l, int _ml, double norm)
: Base(e,l,norm)
, ml(_ml)
{};



}}//namespace
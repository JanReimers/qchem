// File:: BasisSet1/Atom/Evaluators/Imp/Evaluator.C
module;
#include <blaze/Math.h>

module qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.Symmetry.Spherical;

namespace BasisSet::Atom::Evaluators
{
Evaluator::Evaluator(const sym_t& s)
    : l(Symmetry::Getl(s))

    , grouper(0)
{}

std::string Evaluator::AngularID() const
{
    return std::to_string(l);
}

} //namespace
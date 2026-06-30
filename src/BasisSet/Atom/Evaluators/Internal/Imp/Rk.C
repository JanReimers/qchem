// File: BasisSet/Atom/Evaluators/Internal/Imp/Rk.C
module;
#include <cassert>
module qchem.BasisSet.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet.Atom.Evaluators;

namespace qchem {

bool Rk::isSupported(const Cache4_Client* cl) const
{
    auto eval=dynamic_cast<const BasisSet::Atom::Evaluators::Evaluator*>(cl);
    assert(eval);
    return eval->Getl()<=LMax();
}

} // namespace qchem
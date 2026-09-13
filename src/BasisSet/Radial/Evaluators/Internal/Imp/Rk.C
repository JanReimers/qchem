// File: BasisSet/Radial/Evaluators/Internal/Imp/Rk.C
module;
#include <cassert>
module qchem.BasisSet.Radial.Evaluators.Internal.Rk;
import qchem.BasisSet.Radial.Evaluators;

namespace qchem {

bool Rk::isSupported(const Cache4_Client* cl) const
{
    auto eval=dynamic_cast<const BasisSet::Radial::Evaluators::Evaluator*>(cl);
    assert(eval);
    return eval->Getl()<=LMax();
}

} // namespace qchem
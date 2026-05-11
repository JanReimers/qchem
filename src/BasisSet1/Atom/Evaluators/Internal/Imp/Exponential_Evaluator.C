// File: BasisSet1/Atom/Evaluators/Internal/Imp/Exponential_Evaluator.C  Common base for Slater and Gaussian evaluators
module;
#include <cmath>
#include <cassert>
#include <iostream>
#include <blaze/math/DynamicVector.h>

module qchem.BasisSet1.Atom.Evaluators.Internal.Exponential_IBS_Evaluator;

void Exponential_IBS_Evaluator::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<ExponentGrouper*>(_grouper);
    assert(grouper);
    for (auto e:es) es_indices.push_back(grouper->Insert(e,l));
    // std::cout << "es_indices=" << es_indices << std::endl;
}

std::string Exponential_IBS_Evaluator::RadialID () const
{
    std::ostringstream os;
    os << Name() << "{";
    for (auto e:es) os << e << " ";
    os << "}";
    return os.str();
}


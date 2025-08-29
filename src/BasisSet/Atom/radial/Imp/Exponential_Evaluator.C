// File: BasisSet/Atom/radial/Imp/Exponential_Evaluator.C  Common for Slater and Gaussian evaluators
module;
#include <valarray>
#include <cmath>
#include <cassert>
#include <iostream>

module BasisSet.Atom.Exponential_IBS_Evaluator;

void Exponential_IBS_Evaluator::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<ExponentGrouper*>(_grouper);
    assert(grouper);
    for (auto e:es) es_indices.push_back(grouper->Insert(e,l));
    // std::cout << "es_indices=" << es_indices << std::endl;
}

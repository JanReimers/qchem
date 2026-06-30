// File: BasisSet/Atom/Evaluators/Internal/Imp/Exponential_Evaluator.C  Common base for Slater and Gaussian evaluators
module;
#include <cassert>
#include <iostream>
module qchem.BasisSet.Atom.Evaluators.Internal.ExponentialEvaluator;
import qchem.Math;
import qchem.Blaze;

namespace qchem::BasisSet::Atom::Evaluators
{

void ExponentialEvaluator::Register(Grouper* _grouper)
{
    assert(_grouper);
    auto grouper=static_cast<ExponentGrouper*>(_grouper);
    assert(grouper);
    for (auto e:es) es_indices.push_back(grouper->Insert(e,Getl()));
    // std::cout << "es_indices=" << es_indices << std::endl;
}

std::string ExponentialEvaluator::RadialID () const
{
    std::ostringstream os;
    if (isEvenTempered)
    {
        os << Name() << " N=" << es.size() << " {";
        if (es.size()>0) os << es[0];
        if (es.size()>1) os << " ... " << es[es.size()-1];
    }
    else
    {
        os << Name() << " {";
        for (auto e:es) os << e << " ";
    }
    os << "}";
    return os.str();
}

bool ExponentialEvaluator::EvenTempered(const rvec_t& es)
{
    bool et=true;
    if (es.size()>1)
    { 
        double beta=es[1]/es[0];
        for (size_t i=2;i<es.size();i++)
            if (fabs(es[i]/(beta*es[i-1])-1.0)>1e-14)
            {
                std::cout << "Warning: Irrep basis set is not even tempered fabs(es[" << i << "]/(beta*es[" << i-1 << "])-1.0) = " << fabs(es[i]/(beta*es[i-1])-1.0) << std::endl;
                et=false;
                break;
            }
    }
    return et;
}

} //namespace
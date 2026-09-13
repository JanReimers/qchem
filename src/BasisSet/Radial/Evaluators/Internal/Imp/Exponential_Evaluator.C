// File: BasisSet/Radial/Evaluators/Internal/Imp/Exponential_Evaluator.C  Common base for Slater and Gaussian evaluators
module;
#include <cassert>
#include <iostream>
#include <sstream>   // the spatial irrep label for Announce
#include <nlohmann/json.hpp>
module qchem.BasisSet.Radial.Evaluators.Internal.ExponentialEvaluator;
import qchem.Math;
import qchem.Blaze;
import qchem.Reporting;   // serialize the exponents into the run report (a report-only sink, never a getter)

namespace qchem::BasisSet::Radial::Evaluators
{

// Announce this shell's exponents as a basis.exponents row {irrep, values}.  The irrep label is the SPATIAL
// one ("s"/"p"/...): exponents know nothing of spin, so a polarized consumer (CLIapps/valgen) strips the
// spin suffix off its usage rows to join.  Silent unless a "basis" section is open -- the orchestrator's
// choice of context (the same gate the LASolver's conditioning write uses), so a basis built outside any
// run pays nothing and a run that did not build its basis truthfully carries no exponents rows.  The
// exponents are SERIALIZED here, never returned for computation -- they stay encapsulated.
void ExponentialEvaluator::Announce(const sym_t& ir) const
{
    namespace rpt = qchem::report;
    if (!rpt::InSection("basis")) return;
    rpt::json values = rpt::json::array();
    for (auto e : es) values.push_back(double(e));
    std::ostringstream os; os << *ir;
    rpt::Row r("exponents");
    rpt::Set("irrep",  os.str());
    rpt::Set("values", values);
}

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
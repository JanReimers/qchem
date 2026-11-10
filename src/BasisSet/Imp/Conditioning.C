// File: BasisSet/Imp/Conditioning.C  The basis conditioning report.
module;
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <string>
#include <type_traits>
module qchem.BasisSet.Conditioning;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.LASolver;
import qchem.LASolver.OverlapConditioning;
import qchem.Reporting;
import qchem.Symmetry;

namespace qchem::BasisSet
{

//! Every block of \a bs at its CONSTRUCTED scalar (doc/RealComplexPlan.md Step 3c-2): the cross-scalar real view first, else the
//! same-scalar block, so no caller downcasts to the concrete container.
template <class T, class F> static void ForEachBlock(const tBasisSet<T>& bs, F&& f)
{
    for (size_t i=0;i<bs.GetNumIBS();++i)
        if (const auto* r=bs.GetRealIBS(i)) f(r);
        else                                f(bs[i]);
}

template <class T> size_t ReportConditioning(const tBasisSet<T>& bs)
{
    namespace rpt = qchem::report;
    size_t removed=0, wouldDrop=0, nBlocks=0;
    double lambdaMin=0.0, cond=0.0;
    std::string lambdaMinAt, condAt;
    ForEachBlock(bs, [&](const auto* b)
    {
        const auto& S=b->Overlap();   // analytic: no grids needed
        using U=typename std::decay_t<decltype(S)>::ElementType;
        const std::string label=b->GetSymmetry().GetLabel();
        const OverlapConditioning c=AnalyzeOverlap<U>(S);
        {
            rpt::Row row("perIrrep");
            rpt::Set("irrep",      label);
            rpt::Set("nFunctions", (long)b->GetNumFunctions());
            rpt::Set("real",       b->IsReal());   // TRIM fact from the irrep (doc/RealComplexPlan.md Step 1)
            rpt::Set("runsReal",   std::is_same_v<U,double>);   // the block's CONSTRUCTED scalar (3c-3 evidence)
            rpt::Set("lambdaMin",  c.lambdaMin);
            rpt::Set("lambdaMax",  c.lambdaMax);
            rpt::Set("cond",       c.cond);
        }
        for (size_t idx : PivotedCholeskyDrops<U>(S))   // redundant AOs (empty on a healthy basis)
        {
            rpt::Row r("removed");
            rpt::Set("irrep", label);
            rpt::Set("index", (long)idx);
            ++removed;
        }
        if (S.rows()==0) return;
        if (nBlocks==0 || c.lambdaMin<lambdaMin) { lambdaMin=c.lambdaMin; lambdaMinAt=label; }
        if (nBlocks==0 || c.cond>cond)           { cond=c.cond;           condAt=label; }
        wouldDrop+=c.wouldDrop;
        ++nBlocks;
    });
    {
        rpt::Scope summary("conditioning");
        rpt::Set("lambdaMin",   lambdaMin);
        rpt::Set("lambdaMinAt", lambdaMinAt);
        rpt::Set("cond",        cond);
        rpt::Set("condAt",      condAt);
        rpt::Set("floor",       ConditioningFloor);
        rpt::Set("wouldDrop",   (long)wouldDrop);
    }
    // Straight to cout, like the other tagged lines: report::Log is silent in a run with no console, and this line must never be.
    report::FormatGuard keep(std::cout);
    std::cout<<"[basis trim] report only: "<<nBlocks<<" block(s), min eig S "<<std::scientific<<std::setprecision(3)<<lambdaMin<<" at "
             <<lambdaMinAt<<", max cond(S) "<<cond<<" at "<<condAt<<", "<<wouldDrop<<" AO(s) would drop at floor "<<std::setprecision(0)
        <<ConditioningFloor<<std::defaultfloat<<std::endl;
    return removed;
}

template <class T> double MinOverlapEigenvalue(const tBasisSet<T>& bs)
{
    double lmin=1e300;
    ForEachBlock(bs, [&](const auto* b)
    {
        const auto& S=b->Overlap();
        using U=typename std::decay_t<decltype(S)>::ElementType;
        lmin=std::min(lmin, AnalyzeOverlap<U>(S).lambdaMin);
    });
    return lmin;
}

template size_t ReportConditioning<double>(const tBasisSet<double>&);
template size_t ReportConditioning<dcmplx>(const tBasisSet<dcmplx>&);
template double MinOverlapEigenvalue<double>(const tBasisSet<double>&);
template double MinOverlapEigenvalue<dcmplx>(const tBasisSet<dcmplx>&);

} // namespace qchem::BasisSet

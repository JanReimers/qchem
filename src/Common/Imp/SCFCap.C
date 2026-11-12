// File: Imp/SCFCap.C  The cap-hit hook and the bounded-SCF scope (see qchem.SCFCap).
module;
#include <cassert>
#include <cstddef>
#include <functional>
#include <iostream>
#include <string>
module qchem.SCFCap;

namespace qchem
{

namespace {
SCFCapHook& Hook() {static SCFCapHook h; return h;}
thread_local int g_bounded=0;   //!< live BoundedSCF scopes on this thread (the SCF and its test share the thread)
thread_local IterationBudget* g_budget=nullptr;   //!< the innermost live IterationBudget on this thread
thread_local int g_audit=0;   //!< live AuditSCF scopes on this thread
} // namespace

void SetSCFCapHook(SCFCapHook hook) {Hook()=std::move(hook);}

void ReportSCFCapHit(const std::string& what)
{
    if (g_bounded>0 || !Hook()) return;
    Hook()(what);
}

BoundedSCF::BoundedSCF(const char* why)
{
    assert(why && *why && "BoundedSCF: say why this SCF is deliberately bounded");
    (void)why;
    ++g_bounded;
}
BoundedSCF::~BoundedSCF() {--g_bounded;}

void ReportSCFIterations(std::size_t n)
{
    for (IterationBudget* b=g_budget; b; b=b->itsOuter) (g_audit>0 ? b->itsAudit : b->itsUsed)+=n;
}

AuditSCF:: AuditSCF() {++g_audit;}
AuditSCF::~AuditSCF() {--g_audit;}

SetupSCF:: SetupSCF() : itsSaved(g_budget) {g_budget=nullptr;}
SetupSCF::~SetupSCF() {g_budget=itsSaved;}

IterationBudget::IterationBudget(std::size_t budget, const char* what)
    : itsBudget(budget), itsWhat(what), itsOuter(g_budget)
{
    assert(what && *what && "IterationBudget: name the panel case");
    g_budget=this;
}
IterationBudget::~IterationBudget()
{
    assert(g_budget==this && "IterationBudget: scopes must nest");
    g_budget=itsOuter;
    std::cout << "[panel] " << itsWhat << ": " << itsUsed << " SCF iterations (budget " << itsBudget << ")";
    if (itsAudit>0) std::cout << " + " << itsAudit << " audit re-convergence iterations (not budgeted)";
    std::cout << std::endl;
    if (itsUsed>itsBudget && Hook())
        Hook()("[panel] "+std::string(itsWhat)+": "+std::to_string(itsUsed)+" SCF iterations, OVER its budget of "
               +std::to_string(itsBudget)+" (doc/SCFRobustnessPlan.md §4)");
}

} // namespace

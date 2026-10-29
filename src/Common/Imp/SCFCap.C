// File: Imp/SCFCap.C  The cap-hit hook and the bounded-SCF scope (see qchem.SCFCap).
module;
#include <cassert>
#include <functional>
#include <string>
module qchem.SCFCap;

namespace qchem
{

namespace {
SCFCapHook& Hook() {static SCFCapHook h; return h;}
thread_local int g_bounded=0;   //!< live BoundedSCF scopes on this thread (the SCF and its test share the thread)
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

} // namespace

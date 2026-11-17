// File: BasisSet/Imp/RealTableRelease.C  The release phase's scope: it writes what the phase did into the caller's report context.
module;
#include <cstddef>
#include <iostream>
module qchem.BasisSet.RealTableRelease;
import qchem.Reporting;

namespace qchem::BasisSet
{

RealTablesReleased::RealTablesReleased(const RealTableRelease& fit) : itsFit(fit), itsBytes(fit.BeginRealTableRelease())
{
    std::cout << "[response memory] released the ground state's real-block XC tables for the response: "
              << double(itsBytes)/1e9 << " GB" << std::endl;
    qchem::report::Set("realTableBytesReleased", itsBytes);
}

RealTablesReleased::~RealTablesReleased()
{
    const size_t n=itsFit.EndRealTableRelease();
    std::cout << "[response memory] release phase closed: " << n << " real-block table(s) rebuilt inside it" << std::endl;
    qchem::report::Set("realTablesRebuiltInPhase", n);
}

size_t RealTablesReleased::BytesReleased() const {return itsBytes;}

} // namespace

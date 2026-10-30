// File: Common/HandOff.C  THE SCF HAND-OFF LINE (doc/SCFRobustnessPlan.md R1.3, 2026-10-30).
//
// Every point where one extrapolator, mixer regime or stage gives way to another prints ONE line in ONE format, so a run's
// hand-offs can be read (and grepped) off the console without knowing which object made them.  The object that DECIDES a
// hand-off is often not the one that can ACT on it: a ladder rung change is decided inside the accelerator library, which cannot
// see the density mixer.  So a decider REPORTS a HandOffRequest through its own interface, and the owner of both (the iterator,
// or the calculation for a stage) acts and prints.  This module is the shared vocabulary; it lives in qcCommon because both
// qcSCFAccelerator and qcChargeDensity report through it and neither may see the other.
module;
#include <cstddef>
#include <iostream>
#include <string>
export module qchem.HandOff;

export namespace qchem
{
//! One hand-off: what gave way, to what, and why (the decider's own trigger, in words).
struct HandOffRequest
{
    std::string from, to, why;
};

//! THE one line: `[handoff] iteration N: <from> -> <to> (<why>)<note>`.  Iteration 0 = before the stage's first iteration.
inline void PrintHandOff(std::size_t iteration, const HandOffRequest& h, const std::string& note={})
{
    std::cout << "[handoff] iteration " << iteration << ": " << h.from << " -> " << h.to << " (" << h.why << ")" << note << std::endl;
}
} // namespace qchem

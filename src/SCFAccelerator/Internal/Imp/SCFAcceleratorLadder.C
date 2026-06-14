// File: SCFAcceleratorLadder.C  Implementation of the chained SCF accelerator.
module;
#include <iostream>
#include <vector>
module qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;

namespace qchem::SCFAccelerators
{
using std::cout;
using std::endl;

SCFIrrepAcceleratorLadder::~SCFIrrepAcceleratorLadder() { for (auto r:itsRungs) delete r; }

void SCFIrrepAcceleratorLadder::UseFD(const smat_t<double>& F, const smat_t<double>& DPrime)
{
    // Feed every rung so whichever becomes active next already has the current F.
    for (auto r:itsRungs) r->UseFD(F,DPrime);
}
LASolver<double>::UUd_t SCFIrrepAcceleratorLadder::NextOrbitals()
{
    return itsRungs[*itsActive]->NextOrbitals();
}

//----------------------------------------------------------------------------------------
SCFAcceleratorLadder::~SCFAcceleratorLadder() { for (auto r:itsRungs) delete r; }

SCFIrrepAccelerator* SCFAcceleratorLadder::Create(const LASolver<double>* las,const Irrep& qns,int occ)
{
    std::vector<SCFIrrepAccelerator*> rungs;
    for (auto r:itsRungs) rungs.push_back(r->Create(las,qns,occ));
    return new SCFIrrepAcceleratorLadder(std::move(rungs), &itsActive);
}

bool SCFAcceleratorLadder::CalculateProjections()
{
    bool ok = itsRungs[itsActive]->CalculateProjections();

    // Track the best error so far on this rung.  A rung that keeps beating its own best is
    // still making progress (however slowly) -> do NOT hand off.  Only switch when the rung
    // is Exhausted() AND has failed to improve for several consecutive steps (genuinely
    // stuck/oscillating, not merely slow -- a heavy atom can converge slowly yet steadily).
    double err = itsRungs[itsActive]->GetError();
    if (err < 0.999*itsBestErr) { itsBestErr=err; itsNoImprove=0; }
    else                          itsNoImprove++;

    // Error floor: a plateau at the noise floor is convergence, not a stall -- never hand off
    // there (e.g. a heavy atom whose DIIS bottoms out at [F,D]~1e-6 below the SCF criteria).
    const double floor=1e-4;
    if (itsActive+1<itsRungs.size() && itsRungs[itsActive]->Exhausted()
        && itsNoImprove>=5 && err>floor)
    {
        cout << "  *** SCF accelerator ladder: rung " << itsActive
             << " exhausted -> advancing to rung " << itsActive+1 << " ***" << endl;
        itsActive++;
        itsBestErr=1e300; itsNoImprove=0;
    }
    return ok;
}

double SCFAcceleratorLadder::GetError() const { return itsRungs[itsActive]->GetError(); }
void   SCFAcceleratorLadder::ShowLabels(std::ostream& os)      const { itsRungs[itsActive]->ShowLabels(os); }
void   SCFAcceleratorLadder::ShowConvergence(std::ostream& os) const { itsRungs[itsActive]->ShowConvergence(os); }

} //namespace

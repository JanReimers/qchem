// File: SCFAcceleratorLadder.C  Implementation of the chained SCF accelerator.
module;
#include <iostream>
#include <vector>
#include <deque>
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

bool SCFAcceleratorLadder::Stalled() const
{
    // Error has not dropped by 2x across the recent window -> the active rung is stuck.
    if (itsErr.size()<4) return false;
    return itsErr.back() > 0.5*itsErr.front();
}

bool SCFAcceleratorLadder::CalculateProjections()
{
    bool ok = itsRungs[itsActive]->CalculateProjections();
    itsErr.push_back(itsRungs[itsActive]->GetError());
    if (itsErr.size()>5) itsErr.pop_front();

    if (itsActive+1<itsRungs.size() && itsRungs[itsActive]->Exhausted() && Stalled())
    {
        cout << "  *** SCF accelerator ladder: rung " << itsActive
             << " exhausted -> advancing to rung " << itsActive+1 << " ***" << endl;
        itsActive++;
        itsErr.clear();
    }
    return ok;
}

double SCFAcceleratorLadder::GetError() const { return itsRungs[itsActive]->GetError(); }
void   SCFAcceleratorLadder::ShowLabels(std::ostream& os)      const { itsRungs[itsActive]->ShowLabels(os); }
void   SCFAcceleratorLadder::ShowConvergence(std::ostream& os) const { itsRungs[itsActive]->ShowConvergence(os); }

} //namespace

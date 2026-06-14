// File: SCFAcceleratorLadder.C  Implementation of the chained SCF accelerator.
module;
#include <iostream>
#include <vector>
#include <cmath>
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

void SCFAcceleratorLadder::SetEnergy(double E) { itsPrevE=itsLastE; itsLastE=E; }

bool SCFAcceleratorLadder::CalculateProjections()
{
    bool ok = itsRungs[itsActive]->CalculateProjections();

    // Track the best error so far on this rung: a rung that keeps beating its own best is
    // still making progress, so it is not "stuck".
    double err = itsRungs[itsActive]->GetError();
    if (err < 0.999*itsBestErr) { itsBestErr=err; itsNoImprove=0; }
    else                          itsNoImprove++;

    // Hand off only if the active rung is (a) Exhausted(), (b) stuck on its error for `stall`
    // steps, (c) the ENERGY is still moving (not a converged plateau -- see the design notes
    // at the top of the header), and (d) above an absolute noise floor.  See those notes for
    // why the energy, not [F,D], is the decisive signal.
    double relDE = (itsLastE!=0.0) ? std::fabs((itsLastE-itsPrevE)/itsLastE) : 1.0;
    if (itsActive+1<itsRungs.size()
        && itsRungs[itsActive]->Exhausted()
        && itsNoImprove>=itsStall
        && relDE>itsEThresh
        && err>itsFloor)
    {
        cout << "  *** SCF accelerator ladder: rung " << itsActive
             << " exhausted (|dE/E|=" << relDE << ") -> advancing to rung " << itsActive+1 << " ***" << endl;
        itsActive++;
        itsBestErr=1e300; itsNoImprove=0;
    }
    return ok;
}

double SCFAcceleratorLadder::GetError() const { return itsRungs[itsActive]->GetError(); }
void   SCFAcceleratorLadder::ShowLabels(std::ostream& os)      const { itsRungs[itsActive]->ShowLabels(os); }
void   SCFAcceleratorLadder::ShowConvergence(std::ostream& os) const { itsRungs[itsActive]->ShowConvergence(os); }

} //namespace

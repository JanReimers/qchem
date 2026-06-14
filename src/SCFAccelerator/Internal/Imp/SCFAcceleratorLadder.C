// File: SCFAcceleratorLadder.C  Implementation of the chained SCF accelerator.
module;
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
module qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;

namespace qchem::SCFAccelerators
{
using std::cout;
using std::endl;

// Both ladders own their rungs through unique_ptr, so no explicit destructors are needed.

void SCFIrrepAcceleratorLadder::UseFD(const smat_t<double>& F, const smat_t<double>& DPrime)
{
    // Feed every rung so whichever becomes active next already has the current F.
    for (auto& r:itsRungs) r->UseFD(F,DPrime);
}
LASolver<double>::UUd_t SCFIrrepAcceleratorLadder::NextOrbitals()
{
    return Active()->NextOrbitals();
}

//----------------------------------------------------------------------------------------
SCFIrrepAccelerator* SCFAcceleratorLadder::Create(const LASolver<double>* las,const Irrep& qns,int occ)
{
    std::vector<std::unique_ptr<SCFIrrepAccelerator>> rungs;
    for (auto& r:itsRungs) rungs.emplace_back(r->Create(las,qns,occ));
    return new SCFIrrepAcceleratorLadder(std::move(rungs), &itsActive);
}

void SCFAcceleratorLadder::SetEnergy(double E) { itsPrevE=itsLastE; itsLastE=E; }

bool SCFAcceleratorLadder::CalculateProjections()
{
    bool ok = Active()->CalculateProjections();

    // Track the best error so far on this rung: a rung that keeps beating its own best is
    // still making progress, so it is not "stuck".
    double err = Active()->GetError();
    if (err < 0.999*itsBestErr) { itsBestErr=err; itsNoImprove=0; }
    else                          itsNoImprove++;

    // Hand off only if the active rung is (a) Exhausted(), (b) stuck on its error for `stall`
    // steps, (c) the ENERGY is still moving (not a converged plateau -- see the design notes
    // at the top of the header), and (d) above an absolute noise floor.  See those notes for
    // why the energy, not [F,D], is the decisive signal.
    double relDE = (itsLastE!=0.0) ? std::fabs((itsLastE-itsPrevE)/itsLastE) : 1.0;
    bool stallSwitch = Active()->Exhausted()
                    && itsNoImprove>=itsStall && relDE>itsEThresh && err>itsFloor;
    //near convergence: hand to the polisher (err>0 guards the initial state, before the active
    //rung has produced orbitals and computed a real [F,D]).
    bool tailSwitch  = itsSwitchAt>0.0 && err>0.0 && err<itsSwitchAt;
    if (itsActive+1<itsRungs.size() && (stallSwitch || tailSwitch))
    {
        cout << "  *** SCF accelerator ladder: rung " << itsActive
             << (tailSwitch ? " near convergence (err=" : " exhausted (|dE/E|=")
             << (tailSwitch ? err : relDE)
             << ") -> advancing to rung " << itsActive+1 << " ***" << endl;
        itsActive++;
        itsBestErr=1e300; itsNoImprove=0;
    }
    return ok;
}

// The ladder runs the direct-min loop exactly when its active rung is a direct minimizer.
bool SCFAcceleratorLadder::WantsLineSearch() const { return Active()->WantsLineSearch(); }

double SCFAcceleratorLadder::GetError() const { return Active()->GetError(); }
void   SCFAcceleratorLadder::ShowLabels(std::ostream& os)      const { Active()->ShowLabels(os); }
void   SCFAcceleratorLadder::ShowConvergence(std::ostream& os) const { Active()->ShowConvergence(os); }

} //namespace

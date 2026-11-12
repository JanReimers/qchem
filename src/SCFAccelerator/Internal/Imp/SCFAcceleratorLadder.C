// File: SCFAcceleratorLadder.C  Implementation of the chained SCF accelerator (templated on T).
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
module qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;
import qchem.Reporting;   // the accelerator files what it was BUILT with, one level into its builder's context

namespace qchem::SCFAccelerators
{
using std::cout;
using std::endl;

// The rungs filed themselves as they were built, each into a `rungs` row its builder opened; the ladder files its own knobs.
SCFAcceleratorLadder::SCFAcceleratorLadder(std::vector<std::unique_ptr<SCFAccelerator>> rungs, double ethresh, int stall,
                                           double floor, double switchat, ScheduleSignal signal, size_t start)
    : itsRungs(std::move(rungs)), itsEThresh(ethresh), itsStall(stall),
      itsFloor(floor), itsSwitchAt(switchat), itsSignal(signal), itsActive(start)
{
    assert(itsActive<itsRungs.size());
    report::Set("type",      "Ladder");
    report::Set("eThresh",   ethresh);
    report::Set("stall",     stall);
    report::Set("floor",     floor);
    report::Set("switchAt",  switchat);
    report::Set("signal",    signal==ScheduleSignal::EnergyChange ? "EnergyChange" : "Error");
    report::Set("startRung", (long)start);
}

// Both ladders own their rungs through unique_ptr, so no explicit destructors are needed.

template <class T> void tSCFIrrepAcceleratorLadder<T>::UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)
{
    // Feed every rung so whichever becomes active next already has the current F.
    for (auto& r:itsRungs) r->UseFD(F,DPrime);
}
template <class T> typename LASolver<T>::UUd_t tSCFIrrepAcceleratorLadder<T>::NextOrbitals()
{
    return Active()->NextOrbitals();
}

//----------------------------------------------------------------------------------------
// ONE templated body behind the two typed Create overloads (§6): every rung is itself scalar-agnostic,
// so the per-irrep ladder simply collects each rung's block-typed child.
template <class T> tSCFIrrepAccelerator<T>* SCFAcceleratorLadder::CreateT(const LASolver<T>* las,const Irrep& qns,int occ)
{
    std::vector<std::unique_ptr<tSCFIrrepAccelerator<T>>> rungs;
    for (auto& r:itsRungs) rungs.emplace_back(r->Create(las,qns,occ));
    return new tSCFIrrepAcceleratorLadder<T>(std::move(rungs), &itsActive);
}
tSCFIrrepAccelerator<double>* SCFAcceleratorLadder::Create(const LASolver<double>* las,const Irrep& qns,int occ)
{ return CreateT(las,qns,occ); }
tSCFIrrepAccelerator<dcmplx>* SCFAcceleratorLadder::Create(const LASolver<dcmplx>* las,const Irrep& qns,int occ)
{ return CreateT(las,qns,occ); }

void SCFAcceleratorLadder::SetEnergy(double E) { itsPrevE=itsLastE; itsLastE=E; }

bool SCFAcceleratorLadder::CalculateProjections()
{
    bool ok = Active()->CalculateProjections();

    // Track the best error so far on this rung: a rung that keeps beating its own best is
    // still making progress, so it is not "stuck".  The tracker starts at the FIRST REAL error: before
    // any orbitals exist the rung reports err = 0, and a best of 0 is never beaten, so the count used to
    // run from iteration 1 (lane A, 2026-11-10 -- masked while the stall also required Exhausted()).
    double err = Active()->GetError();
    if      (err<=0.0)                 {}
    else if (err < 0.999*itsBestErr) { itsBestErr=err; itsNoImprove=0; }
    else                               itsNoImprove++;

    // STALL: hand off if the active rung's error has not beaten its own best for `stall` steps, while
    // (c) the ENERGY is still moving (not a converged plateau -- see the design notes at the top of the
    // header) and (d) the error is above an absolute noise floor.  NOT gated on Exhausted() (user ruling
    // 2026-11-10, "yes that sounds like a good exit policy"): a DIIS LIMIT CYCLE extrapolates fine every step
    // and goes nowhere, so it is never Exhausted (B singular), and the tail signal fires on it only by luck.
    // MEASURED on GPW_Si.Γ_U_Ladder_LinearResponseUIsAFixedPoint_Panel (Si +U, DIIS period ~20, [F,D] floor
    // 1.3e-3): the |ΔE/E| tail fired on a lucky tie at iteration 66 and GDM met the 100 cap; this exit hands
    // off at iteration 8 and the run converges in 68 (all three SCFs 90 vs 136, the same U fixed point).  On
    // 12 other Ladder cases it changed nothing or saved one iteration.  The alternatives measured there --
    // the tail on [F,D] < 1e-2 / 1e-3, or on both signals -- either fire on a lucky trough or never fire,
    // and cost NaF 5-17 iterations ("Ladder |ΔE/E| hand-off" row, doc/Records/OpenWork_History5.md).
    double relDE = (itsLastE!=0.0) ? std::fabs((itsLastE-itsPrevE)/itsLastE) : 1.0;
    bool stallSwitch = itsNoImprove>=itsStall && relDE>itsEThresh && err>itsFloor;
    // TAIL hand-off: the SCHEDULE SIGNAL has dropped below switchat (near convergence).  Error ([F,D],
    // molecular default) or EnergyChange (|ΔE/E|, solids -- see ScheduleSignal).  The validity guard skips
    // the initial state: Error needs a real [F,D] (err>0, before any orbitals err==0); EnergyChange needs
    // two energies (itsPrevE set), else relDE's sentinel 1.0 would look "not yet settled" anyway.
    const bool energyGated = (itsSignal==ScheduleSignal::EnergyChange);
    const double tailMetric = energyGated ? relDE : err;
    const bool   tailValid  = energyGated ? (itsPrevE!=0.0) : (err>0.0);
    bool tailSwitch  = itsSwitchAt>0.0 && tailValid && tailMetric<itsSwitchAt;
    if (itsActive+1<itsRungs.size() && (stallSwitch || tailSwitch))
    {
        // ENGAGEABILITY VETO: never advance onto a rung whose STANDING precondition the run's occupation
        // structure fails (a Fermi-smeared D' is outside GDM's integer-occupation manifold -- knowable from
        // UseFD, which feeds every rung each iteration).  Without this, the loop degrades to bare mixed
        // fixed-point steps under the new rung's tag, and on an ionic cell that is the unstable iteration
        // the CURRENT rung was taming (MnO run 31: -61.414 -> a -60.49 limit cycle).  Announced once.
        if (!itsRungs[itsActive+1]->Engageable())
        {
            if (!itsVetoAnnounced)
                cout << "  *** SCF accelerator ladder: rung " << itsActive+1 << " ("
                     << itsRungs[itsActive+1]->Tag() << ") is NOT ENGAGEABLE for this run's occupation "
                        "structure (smeared/MOM density outside its manifold) -- staying on rung "
                     << itsActive << " (" << Active()->Tag() << ") ***" << endl;
            itsVetoAnnounced=true;
        }
        else
        {
            std::ostringstream why;
            why << "ladder rung " << itsActive
                << (tailSwitch ? (energyGated ? " near convergence, |dE/E|=" : " near convergence, err=")
                               : " stalled: no better error in "+std::to_string(itsNoImprove)+" steps, |dE/E|=")
                << (tailSwitch ? tailMetric : relDE);
            itsHandOff=HandOffRequest{Active()->Tag(), itsRungs[itsActive+1]->Tag(), why.str()};
            itsActive++;
            itsBestErr=1e300; itsNoImprove=0; itsVetoAnnounced=false;
        }
    }
    // RETREAT: the ACTIVE rung's standing precondition failed (discovered late -- e.g. GDM's leading-block
    // check needs its own orbitals, so a MOM mismatch only surfaces after the rung has seeded).  Sitting on
    // an unengageable direct-min rung means mixer-only steps forever; fall back to the previous rung, loudly.
    if (itsActive>0 && !Active()->Engageable())
    {
        itsHandOff=HandOffRequest{Active()->Tag(), itsRungs[itsActive-1]->Tag(),
                                  "ladder RETREAT: rung "+std::to_string(itsActive)+" can no longer engage (occupation structure outside its manifold)"};
        itsActive--;
        itsBestErr=1e300; itsNoImprove=0; itsVetoAnnounced=false;
    }
    return ok;
}

// The ladder runs the direct-min loop exactly when its active rung is a direct minimizer.
bool SCFAcceleratorLadder::WantsLineSearch() const { return Active()->WantsLineSearch(); }
bool SCFAcceleratorLadder::CanLineSearch()  const { return Active()->CanLineSearch(); }
bool SCFAcceleratorLadder::RejectStep()        { return Active()->RejectStep(); }

double SCFAcceleratorLadder::GetError() const { return Active()->GetError(); }
void   SCFAcceleratorLadder::ShowLabels(std::ostream& os)      const { Active()->ShowLabels(os); }
void   SCFAcceleratorLadder::ShowConvergence(std::ostream& os) const { Active()->ShowConvergence(os); }
const char* SCFAcceleratorLadder::Tag()   const { return Active()->Tag(); }   // the live rung
int    SCFAcceleratorLadder::Count() const { return Active()->Count(); }
// The rung live NOW is the one whose NextOrbitals made this iterate: a hand-off or retreat happens inside CalculateProjections,
// BEFORE the per-irrep NextOrbitals, which already delegate to the new rung.
bool   SCFAcceleratorLadder::Extrapolated() const { return Active()->Extrapolated(); }
std::vector<std::pair<size_t,double>> SCFAcceleratorLadder::ExtrapolationWeights() const { return Active()->ExtrapolationWeights(); }
double SCFAcceleratorLadder::MinSV() const { return Active()->MinSV(); }

template class tSCFIrrepAcceleratorLadder<double>;
template class tSCFIrrepAcceleratorLadder<dcmplx>;

} //namespace

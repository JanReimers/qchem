// File: SCFAcceleratorLadder.C  Implementation of the chained SCF accelerator (templated on T).
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
            cout << "  *** SCF accelerator ladder: rung " << itsActive
                 << (tailSwitch ? (energyGated ? " near convergence (|dE/E|=" : " near convergence (err=")
                                : " exhausted (|dE/E|=")
                 << (tailSwitch ? tailMetric : relDE)
                 << ") -> advancing to rung " << itsActive+1 << " ***" << endl;
            itsActive++;
            itsBestErr=1e300; itsNoImprove=0; itsVetoAnnounced=false;
        }
    }
    // RETREAT: the ACTIVE rung's standing precondition failed (discovered late -- e.g. GDM's leading-block
    // check needs its own orbitals, so a MOM mismatch only surfaces after the rung has seeded).  Sitting on
    // an unengageable direct-min rung means mixer-only steps forever; fall back to the previous rung, loudly.
    if (itsActive>0 && !Active()->Engageable())
    {
        cout << "  *** SCF accelerator ladder: rung " << itsActive << " (" << Active()->Tag()
             << ") can no longer engage (occupation structure outside its manifold) -- RETREATING to rung "
             << itsActive-1 << " (" << itsRungs[itsActive-1]->Tag() << ") ***" << endl;
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
double SCFAcceleratorLadder::MinSV() const { return Active()->MinSV(); }

template class tSCFIrrepAcceleratorLadder<double>;
template class tSCFIrrepAcceleratorLadder<dcmplx>;

} //namespace

// File: SCFAcceleratorLadder.C  Chain of SCF accelerators (e.g. DIIS -> GDM).
//
// A composite accelerator: it IS-A SCFAccelerator, so the WaveFunction/SCFIterator see one
// accelerator and never know there is a queue behind it.  It runs the active rung until that
// rung hands off, then advances to the next rung -- which re-seeds itself (e.g. GDM
// diagonalizes on its first step).  Hand-off is hot because the NextOrbitals() interface
// makes every rung an interchangeable orbital producer.
//
//==========================================================================================
// WHEN TO HAND OFF -- design notes (atoms; revisit for molecules/solids/post-HF)
//==========================================================================================
// Deciding *when* a rung is "out of steam" turned out to be the whole problem.  What we
// learned tuning DIIS->GDM on hard cases (heavy DHF atoms, lanthanide BSpline-High HF):
//
//  1. "[F,D] stopped dropping fast" is NOT a usable hand-off signal.  Near convergence DIIS
//     prunes its subspace (singular B, "B.rows()<=2") and [F,D] plateaus at a noise floor
//     while *still converged* -- e.g. Xe DHF bottoms at [F,D]~1e-6, a lanthanide HF stalls
//     ~1e-5.  Reading that plateau as a stall and switching to GDM wrecked the converged Xe
//     result.  An ABSOLUTE [F,D] floor can't separate the two regimes: they are an order of
//     magnitude apart and both oscillate.
//
//  2. The robust signal is the ENERGY, not [F,D].  At every one of those plateaus the energy
//     was already settled (dE/E ~ 1e-13) -- chasing [F,D] further is chasing numerical dust,
//     so GDM there is marginal value and real risk.  GDM's genuine payoff is when DIIS is
//     stuck *and the energy is still moving* (genuine non-convergence: bad guess, near
//     degeneracy, open shell).  So we gate the hand-off on |dE/E| > EThresh AND the rung
//     being Exhausted() AND its error not improving for `stall` steps.
//
//  3. Even when GDM legitimately takes over, its speed-up is currently capped because GDM
//     runs as a per-macro-iteration orbital producer INSIDE the SCF density-mixing loop:
//     the outer relax/MixIn damping fights GDM's direct orbital step and it oscillates in
//     the tail.  The full GDM payoff needs GDM to own the convergence (no outer mixing) with
//     a real geodesic line search -- a separate piece of work.
//
// For molecules/solids/post-HF the *signals* (energy stability + a genuine stall) should
// carry over, but the thresholds (EThresh, stall) will likely need to be revisited per
// problem class, which is why they are run-time params (factory JSON), not constants.
//==========================================================================================
module;
#include <iosfwd>
#include <vector>
export module qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;
export import qchem.SCFAccelerator;

export namespace qchem::SCFAccelerators
{

// Per-irrep: one rung accelerator per ladder rung; delegates to the active one.
class SCFIrrepAcceleratorLadder : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAcceleratorLadder(std::vector<SCFIrrepAccelerator*> rungs, const size_t* active)
        : itsRungs(std::move(rungs)), itsActive(active) {}
    virtual ~SCFIrrepAcceleratorLadder();
    virtual void UseFD(const smat_t<double>& F, const smat_t<double>& DPrime);
    virtual LASolver<double>::UUd_t NextOrbitals();
private:
    std::vector<SCFIrrepAccelerator*> itsRungs;
    const size_t*                     itsActive; //shared with the top-level ladder
};

// Top-level: chain {DIIS, GDM, ...}.  Switch when the active rung is Exhausted() and stalled.
class SCFAcceleratorLadder : public virtual SCFAccelerator
{
public:
    // Hand off when the active rung is Exhausted() AND its error has not improved for `stall`
    // steps AND the energy is still moving (|dE/E| > ethresh).  `floor` is a low backstop:
    // never hand off once the error is below it (an absolute noise floor).
    SCFAcceleratorLadder(std::vector<SCFAccelerator*> rungs,
                         double ethresh=1e-8, int stall=5, double floor=1e-8)
        : itsRungs(std::move(rungs)), itsEThresh(ethresh), itsStall(stall), itsFloor(floor) {}
    virtual ~SCFAcceleratorLadder();
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
    virtual void   SetEnergy(double E);
private:
    std::vector<SCFAccelerator*> itsRungs;
    double                       itsEThresh; //hand off only while |dE/E| exceeds this
    int                          itsStall;
    double                       itsFloor;
    size_t                       itsActive=0;
    double                       itsBestErr=1e300; //best (smallest) error since this rung started
    int                          itsNoImprove=0;   //consecutive steps without beating itsBestErr
    double                       itsLastE=0.0, itsPrevE=0.0; //last two reported total energies
};

} //namespace

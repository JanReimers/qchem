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
//  3. GDM only pays off when it OWNS the loop.  Running GDM as a per-macro-iteration orbital
//     producer inside the SCF density-mixing loop caps it: the outer relax/MixIn damping
//     fights its direct orbital step and it oscillates in the tail.  The fix (now in place):
//     when the active rung WantsLineSearch(), the SCFIterator switches to its direct-min loop
//     -- a geodesic line search, NO density mixing.  Combined with a TAIL hand-off (switchat:
//     once DIIS drives [F,D] below a small threshold, near convergence, hand to the GDM
//     polisher) this is the production recipe.  Measured on Nd (Z=60) BSpline-High HF to
//     [F,D]<1e-8: DIIS-alone takes 57 iters with an oscillating tail (bounces 1e-8<->4e-8);
//     DIIS->direct-min (switchat=1e-4) reaches the SAME energy in 27 iters with a smooth,
//     strictly monotone tail.  (GDM far from convergence is still useless -- direct-min from
//     a cold start on an open shell crawls and stalls -- which is exactly why it is wired as
//     the near-convergence polisher, not a from-scratch solver.)
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
    // Delegate the direct-min line-search hooks to the active rung (so a GDM rung's geodesic
    // step is reachable once the ladder hands off to it).
    virtual bool ComputeStep() { return itsRungs[*itsActive]->ComputeStep(); }
    virtual LASolver<double>::UUd_t OrbitalsAt(double t, bool commit)
        { return itsRungs[*itsActive]->OrbitalsAt(t,commit); }
private:
    std::vector<SCFIrrepAccelerator*> itsRungs;
    const size_t*                     itsActive; //shared with the top-level ladder
};

// Top-level: chain {DIIS, GDM, ...}.  Switch when the active rung is Exhausted() and stalled.
class SCFAcceleratorLadder : public virtual SCFAccelerator
{
public:
    // Two distinct hand-off triggers (see the design notes above):
    //   STALL  -- the active rung is Exhausted() AND its error has not improved for `stall`
    //             steps AND the energy is still moving (|dE/E| > ethresh): a genuine
    //             non-convergence the next rung should rescue.  `floor` is a low backstop.
    //   TAIL   -- the active rung has driven the error below `switchat` (i.e. we are near
    //             convergence): hand off to the next rung to POLISH the tail.  This is the
    //             slot for a direct minimizer (GDM owns the loop), which is fast and robust
    //             near the minimum but useless far from it.  switchat<=0 disables this.
    SCFAcceleratorLadder(std::vector<SCFAccelerator*> rungs,
                         double ethresh=1e-8, int stall=5, double floor=1e-8, double switchat=0.0)
        : itsRungs(std::move(rungs)), itsEThresh(ethresh), itsStall(stall),
          itsFloor(floor), itsSwitchAt(switchat) {}
    virtual ~SCFAcceleratorLadder();
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
    virtual void   SetEnergy(double E);
    virtual bool   WantsLineSearch() const; //true once the active rung is a direct minimizer
private:
    std::vector<SCFAccelerator*> itsRungs;
    double                       itsEThresh; //hand off only while |dE/E| exceeds this
    int                          itsStall;
    double                       itsFloor;
    double                       itsSwitchAt; //tail hand-off: switch once error < this
    size_t                       itsActive=0;
    double                       itsBestErr=1e300; //best (smallest) error since this rung started
    int                          itsNoImprove=0;   //consecutive steps without beating itsBestErr
    double                       itsLastE=0.0, itsPrevE=0.0; //last two reported total energies
};

} //namespace

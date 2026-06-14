// File: SCFAcceleratorLadder.C  Chain of SCF accelerators (e.g. DIIS -> GDM).
//
// A composite accelerator: it IS-A SCFAccelerator, so the WaveFunction/SCFIterator see one
// accelerator and never know there is a queue behind it.  It runs the active rung until that
// rung is Exhausted() AND the error has stalled, then advances to the next rung -- which
// re-seeds itself (e.g. GDM diagonalizes on its first step).  Hand-off is hot because the
// NextOrbitals() interface makes every rung an interchangeable orbital producer.
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
    // floor: never hand off once the error is below this (a noise-floor plateau is convergence).
    // stall: number of no-improvement steps a rung must show before handing off.
    SCFAcceleratorLadder(std::vector<SCFAccelerator*> rungs, double floor=1e-4, int stall=5)
        : itsRungs(std::move(rungs)), itsFloor(floor), itsStall(stall) {}
    virtual ~SCFAcceleratorLadder();
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
private:
    std::vector<SCFAccelerator*> itsRungs;
    double                       itsFloor;
    int                          itsStall;
    size_t                       itsActive=0;
    double                       itsBestErr=1e300; //best (smallest) error since this rung started
    int                          itsNoImprove=0;   //consecutive steps without beating itsBestErr
};

} //namespace

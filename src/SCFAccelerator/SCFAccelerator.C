// FIle: SCFAccelerator.C  Interface for an SCF accelerator alogrithm
module;
#include <iosfwd>
export module qchem.SCFAccelerator;
export import qchem.Symmetry.Irrep;
export import qchem.LASolver;

export namespace qchem::SCFAccelerators
{

class SCFIrrepAccelerator
{
public:
    virtual ~SCFIrrepAccelerator() {};
    // Feed the current (AO) Fock matrix and the (orthonormal-basis) density matrix.
    virtual void UseFD(const smat_t<double>& F, const smat_t<double>& DPrime)=0;
    // Produce the next set of orbital coefficients (U, U', e), as SolveOrtho returns.
    // A Fock-extrapolator (DIIS) extrapolates F then diagonalizes; a direct minimizer
    // (GDM) rotates the current orbitals along the Grassmann manifold instead.
    virtual LASolver<double>::UUd_t NextOrbitals()=0;

    // Direct-minimization line-search hooks (only GDM-style minimizers implement these).
    //   ComputeStep(): compute the search direction/geodesic for the current Fock without
    //     moving the orbitals.  Returns false if the accelerator does not support a line
    //     search, or wants the caller to fall back to NextOrbitals() (e.g. its seed step).
    //   OrbitalsAt(t,commit): orbitals at geodesic fraction t; commit=false is a pure trial
    //     (for evaluating the energy), commit=true takes the step.
    virtual bool ComputeStep() {return false;}
    virtual LASolver<double>::UUd_t OrbitalsAt(double t, bool commit) {return NextOrbitals();}
};

class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep&, int occ)=0;
    virtual bool CalculateProjections()=0;
    virtual void ShowLabels     (std::ostream&) const=0;
    virtual void ShowConvergence(std::ostream&) const=0;
    virtual double GetError() const=0;
    // Has this accelerator run out of steam (ladder hand-off signal)?  Default: never.
    virtual bool Exhausted() const {return false;}
    // The SCF iterator reports the current total energy each macro-iteration.  The ladder
    // uses the energy change to decide hand-offs (see SCFAcceleratorLadder); others ignore it.
    virtual void SetEnergy(double E) {}
    // Should the SCF iterator run its direct-minimization loop (geodesic line search, no
    // density mixing)?  True for a GDM-style minimizer; the ladder reports its active rung's.
    virtual bool WantsLineSearch() const {return false;}
};

} //namespace


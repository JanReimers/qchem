// FIle: SCFAccelerator.C  Interface for an SCF accelerator alogrithm
module;
#include <iosfwd>
export module qchem.SCFAccelerator;
export import qchem.Symmetry.Irrep;
export import qchem.LASolver;

export namespace qchem::SCFAccelerators
{

// Templated on the matrix element type T (rX/cX convention).  For T=double, hmat_t<double> IS
// smat_t<double>=rsmat_t, so the existing real accelerators (DIIS/GDM/Ladder/Null) -- which bind the
// <double> aliases below -- are unchanged.  The plane-wave path uses tSCFAcceleratorNull<dcmplx>.
template <class T> class tSCFIrrepAccelerator
{
public:
    virtual ~tSCFIrrepAccelerator() {};
    // Feed the current (AO) Fock matrix and the (orthonormal-basis) density matrix.
    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime)=0;
    // Produce the next set of orbital coefficients (U, U', e), as SolveOrtho returns.
    // A Fock-extrapolator (DIIS) extrapolates F then diagonalizes; a direct minimizer
    // (GDM) rotates the current orbitals along the Grassmann manifold instead.
    virtual typename LASolver<T>::UUd_t NextOrbitals()=0;

    // Direct-minimization line-search hooks (only GDM-style minimizers implement these).
    //   ComputeStep(): compute the search direction/geodesic for the current Fock without
    //     moving the orbitals.  Returns false if the accelerator does not support a line
    //     search, or wants the caller to fall back to NextOrbitals() (e.g. its seed step).
    //   OrbitalsAt(t,commit): orbitals at geodesic fraction t; commit=false is a pure trial
    //     (for evaluating the energy), commit=true takes the step.
    virtual bool ComputeStep() {return false;}
    virtual typename LASolver<T>::UUd_t OrbitalsAt(double t, bool commit) {return NextOrbitals();}
};

template <class T> class tSCFAccelerator
{
public:
    virtual ~tSCFAccelerator() {};
    virtual tSCFIrrepAccelerator<T>* Create(const LASolver<T>*,const Irrep&, int occ)=0;
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

using SCFIrrepAccelerator  = tSCFIrrepAccelerator<double>;  using cSCFIrrepAccelerator = tSCFIrrepAccelerator<dcmplx>;
using SCFAccelerator       = tSCFAccelerator<double>;       using cSCFAccelerator      = tSCFAccelerator<dcmplx>;

} //namespace


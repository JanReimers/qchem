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
};

} //namespace


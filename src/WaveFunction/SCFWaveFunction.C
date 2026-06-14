// File: SCFWaveFunction.C  The SCF-driver face of a wave function.
//
// Interface Segregation: WaveFunction is the const, queryable "what is the electronic state"
// interface that every client uses (testers, persistence, future property/post-HF code).
// SCFWaveFunction adds the mutating, iteration-driving methods that ONLY the SCFIterator
// calls.  A wave function that is not produced by an SCF loop (e.g. deserialized from disk,
// or a future post-HF result) implements WaveFunction but NOT SCFWaveFunction.
//
// Virtual inheritance: CompositeWF and Un/PolarizedWF form a diamond on WaveFunction, so
// SCFWaveFunction inherits it virtually to keep a single shared WaveFunction subobject.
module;
export module qchem.WaveFunction.SCF;
export import qchem.WaveFunction;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.Orbitals;
export import qchem.Symmetry.Irrep;

namespace qchem::WaveFunction
{

using Hamiltonian::Hamiltonian;
using ChargeDensity::DM_CD;
using Orbitals::Orbitals;

export class SCFWaveFunction : public virtual WaveFunction
{
public:
    using WaveFunction::GetOrbitals; //keep the const overload visible alongside the mutable one

    virtual void       DoSCFIteration  (Hamiltonian&,const DM_CD*)              =0;
    // Direct-minimization hooks (cf. the SCFIterator direct-min loop):
    //   build the Fock and ask each accelerator to compute its step (no orbital move);
    //   returns false in the seed step (the caller should DoSCFIteration to diagonalize).
    virtual bool       BuildFockAndComputeSteps(Hamiltonian&,const DM_CD*)      =0;
    //   move the orbitals to geodesic fraction t (commit=false is a line-search trial) and refill.
    virtual void       MoveOrbitals    (double t, bool commit, double mergeTol) =0;
    virtual void       FillOrbitals    (double mergeTol)                        =0; //WF knows the electronic structure
    virtual Orbitals*  GetOrbitals     (const Irrep&)                           =0; //mutable access for the loop
};

} //namespace

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
using qchem::Hamiltonian::tHamiltonian;  // qchem:: qualifies the namespace (the alias above shadows it)
using ChargeDensity::DM_CD;
using ChargeDensity::tDM_CD;
using ChargeDensity::tChargeDensity;
using Orbitals::Orbitals;

export template <class T> class tSCFWaveFunction : public virtual tWaveFunction<T>
{
public:
    using tWaveFunction<T>::GetOrbitals; //keep the const overload visible alongside the mutable one

    virtual void       DoSCFIteration  (tHamiltonian<T>&,const tChargeDensity<T>*)      =0;
    //! Iteration-0 seed step: build the Fock from a \a seed density, diagonalize, fill, and return the FIRST
    //! real (matrix-backed) density.  Bundles the DoSCFIteration+FillOrbitals+GetChargeDensity the SCFIterator
    //! used to do by hand; the tDM_CD return makes explicit that a numeric/fit seed yields a genuine density
    //! matrix (the seam where a ScalarFunction seed / HF-bootstrap will land -- see project_numericcd_refactor).
    virtual tDM_CD<T>* Init            (tHamiltonian<T>&,const tChargeDensity<T>* seed, double mergeTol) =0;
    // Direct-minimization hooks (cf. the SCFIterator direct-min loop):
    //   build the Fock and ask each accelerator to compute its step (no orbital move);
    //   returns false in the seed step (the caller should DoSCFIteration to diagonalize).
    virtual bool       BuildFockAndComputeSteps(tHamiltonian<T>&,const tChargeDensity<T>*) =0;
    //   move the orbitals to geodesic fraction t (commit=false is a line-search trial) and refill.
    virtual void       MoveOrbitals    (double t, bool commit, double mergeTol) =0;
    virtual void       FillOrbitals    (double mergeTol)                        =0; //WF knows the electronic structure
    virtual Orbitals*  GetOrbitals     (const Irrep&)                           =0; //mutable access for the loop
};

export using SCFWaveFunction  = tSCFWaveFunction<double>;
export using cSCFWaveFunction = tSCFWaveFunction<dcmplx>;

} //namespace

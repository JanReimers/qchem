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
#include <memory>   // unique_ptr: Init BUILDS the first density (V1.25)
#include <vector>
export module qchem.WaveFunction.SCF;
export import qchem.WaveFunction;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.Orbitals;
export import qchem.Symmetry.Irrep;
export import qchem.ElectronConfiguration.OccupationPolicy;  // the fills consume the iterator's policy (V1.11 inc 3)

namespace qchem::WaveFunction
{

using Hamiltonian::tHamiltonian;
using ChargeDensity::rDM_CD;
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
    //! Iteration-0 seed: build the Fock from \a seed, diagonalize, fill, and hand back the first real
    //! density.  BUILDS it, hence the owning return (V1.25 -- it used to be a raw pointer the caller had to
    //! know to adopt).
    //! \a pol is the run's occupation policy (the SCFIterator's slot).  At seed time it is in its DEFAULT
    //! state -- prescribed integer fill, kT=0, no MOM -- which is what a seed fill SHOULD do (no
    //! self-consistent spectrum exists yet to redistribute over); the old implicit version of this fact
    //! lost charge on metals when anything assumed otherwise (the D11 lesson).
    virtual std::unique_ptr<tDM_CD<T>> Init(tHamiltonian<T>&,const tChargeDensity<T>* seed,
                                            OccupationPolicy<T>& pol, double mergeTol)=0;
    // Direct-minimization hooks (cf. the SCFIterator direct-min loop):
    //   build the Fock and ask each accelerator to compute its step (no orbital move);
    //   returns false in the seed step (the caller should DoSCFIteration to diagonalize).
    virtual bool       BuildFockAndComputeSteps(tHamiltonian<T>&,const tChargeDensity<T>*) =0;
    //   move the orbitals to geodesic fraction t (commit=false is a line-search trial) and refill.
    //! Move the orbitals to geodesic fraction \a t (\a commit=false is a line-search trial) and REFILL
    //! under \a pol.  The direct-min driver passes a \c HeldOccupationPolicy (keep the stored occupied
    //! block) for its trials -- a POLICY, not a flag (V1.11 inc 5; the old \c holdBlock bool is gone).
    virtual void       MoveOrbitals    (OccupationPolicy<T>& pol, double t, bool commit, double mergeTol) =0;
    virtual void       FillOrbitals    (OccupationPolicy<T>&, double mergeTol)  =0;
    virtual Orbitals*  GetOrbitals     (const Irrep&)                           =0; //mutable access for the loop

    // NB: the occupation face is GONE from the wave function (V1.11 increment 3; SCFStrategyPlan §5b).
    // SetMOM / SetSmearing / GetEntropyTerm / AdoptMOMReference / ReleaseMOMReference were post-ctor
    // configuration and side-channel state -- every concrete WF defended against the un-configured
    // window, and the seed fill lived inside it (the D11 charge loss).  All of it now lives on the
    // SCFIterator's OccupationPolicy slot: configuration arrives via OccupationPolicy::Configure, the
    // MOM references (grid-continuation adoption, the 0h guard's release) are policy state keyed by
    // Irrep, and the run's Mermin −TS is read from OccupationPolicy::EntropyTerm() after a fill.
};

export using SCFWaveFunction  = tSCFWaveFunction<double>;
export using cSCFWaveFunction = tSCFWaveFunction<dcmplx>;

} //namespace

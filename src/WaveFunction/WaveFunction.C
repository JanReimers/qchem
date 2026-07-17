// File: WaveFunction.C  Interface for a wave function.
module;
#include <vector>
export module qchem.WaveFunction;
export import qchem.EnergyLevel;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;
export import qchem.ElectronConfiguration;
import qchem.ScalarFunction;
export import qchem.Orbitals;


namespace qchem::WaveFunction
{

export using qchem::ChargeDensity::rDM_CD;
export using qchem::ChargeDensity::tDM_CD;
export using qchem::Orbitals::EnergyLevels;
using Orbitals::Orbitals; //Keep this one last, otherwise it interferes with the two previous declarations!

//! Process-wide MOM (Maximum Overlap Method) toggle (default OFF), mirroring the diagnostic flags
//! `ReportGridCharge`/`ReportBandGap`.  When true, the aufbau occupation follows the MAX OVERLAP onto
//! the previous iteration's occupied subspace (captured from iteration 0's seed diagonalization)
//! instead of the eigenvalue order -- occupied-subspace continuity.  This is the fix for the NaF Γ
//! instability (doc/GPWPlan §0b″): a giant-response diffuse virtual periodically dives below the Fermi
//! edge, and plain aufbau captures it (occupation swap → energy spike); MOM keeps the physical
//! {F 2s, F 2p} manifold occupied through the crossing.  Unlike the parked molecular MOM heuristic
//! (which activated only when the DIIS accelerator engaged), this activates as soon as a reference
//! exists, so it works with the NaF Null accelerator.  Flip it around Iterate and reset it.
export bool& EnableMOM() { static bool on = false; return on; }

//! IMOM (Initial-MOM) reference-capture iteration.  MOM needs a GOOD reference -- the physical occupied
//! subspace -- and the seed (iteration 0) is mid-transient (shapes still shifting fast), so anchoring
//! there locks onto garbage.  Instead run PLAIN aufbau for this many fills to let the SCF descend to the
//! physical fixed point, THEN capture the occupied subspace ONCE and hold it fixed for the rest of the run
//! (a diving diffuse virtual then stays out of the occupied set for good).  Default 10 (NaF descends
//! cleanly by ~iter 8-13, before the first occupation-swap spike at ~14).
export int& MOMStartIter() { static int n = 10; return n; }

// The const, queryable interface to a wave function: "what is the electronic state".
// Every client (testers, persistence, future property/post-HF code) depends on this.
// The mutating, SCF-loop-driving methods live in SCFWaveFunction (see SCFWaveFunction.C),
// which only the SCFIterator uses -- an Interface Segregation split.
//
// Templated on the matrix element type T (rX/cX); WaveFunction is the <double> alias (atoms/
// molecules), cWaveFunction the <dcmplx> instantiation (plane-wave / Bloch-irrep) -- only the
// charge-density type varies with T (the spin density stays a real ScalarFunction).
export template <class T> class tWaveFunction
{
public:
    typedef ScalarFunction<double> sf_t;
    typedef std::vector<Irrep> iqns_t;
    virtual ~tWaveFunction() {};

    virtual const Orbitals* GetOrbitals     (const Irrep&         ) const=0;
    virtual tDM_CD<T>*      GetChargeDensity() const=0;
    virtual sf_t*           GetSpinDensity  () const=0; //Returns a null ptr for un polarized WF.
    virtual EnergyLevels    GetEnergyLevels () const=0;
    virtual iqns_t          GetQNs          () const=0;
    virtual void            DisplayEigen    () const=0;


private:
    tWaveFunction& operator=(const tWaveFunction&);
};

export using WaveFunction  = tWaveFunction<double>;
export using cWaveFunction = tWaveFunction<dcmplx>;

} //namespace

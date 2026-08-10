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

    //! Emit the run report's `basis.usage` section: per-function occupation-weighted populations (the basis-
    //! usage heat map for valence-window tuning, doc/GPWPlan1 §1).  Verbose-only on the console; always recorded
    //! in the json.  Default no-op -- only the composite (atom/molecular) WF implements it (a deserialized or
    //! GPW-Bloch WF opts out), and it is a no-op anyway when no run is open.
    virtual void            EmitBasisUsage  () const {}

private:
    tWaveFunction& operator=(const tWaveFunction&);
};

export using WaveFunction  = tWaveFunction<double>;
export using cWaveFunction = tWaveFunction<dcmplx>;

} //namespace

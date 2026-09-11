// File: WaveFunction.C  Interface for a wave function.
module;
#include <memory>   // std::unique_ptr -- GetChargeDensity BUILDS its density (V1.25)
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
    //! \brief BUILDS the whole-system density -- ALLOCATES, hence the owning return (V1.25).
    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity() const=0;
    virtual EnergyLevels    GetEnergyLevels () const=0;
    virtual iqns_t          GetQNs          () const=0;
    virtual void            DisplayEigen    () const=0;
    // (No Emit*() here.  V1.14: a class does not tell another class WHEN to report -- the composite WF
    // announces its basis usage itself, from FillOrbitals, the moment the occupations exist.)

private:
    tWaveFunction& operator=(const tWaveFunction&);
};

export using WaveFunction  = tWaveFunction<double>;
export using cWaveFunction = tWaveFunction<dcmplx>;

//---------------------------------------------------------------------------------------
//
//  Capability face: a COLLINEAR SPIN-POLARIZED wave function -- the one that HAS a magnetization
//  \f$m(r)=\rho_\uparrow-\rho_\downarrow\f$.  This used to be a pure virtual on tWaveFunction
//  returning a null pointer for the unpolarized half of the hierarchy (V1.17), which declared a
//  capability that half the implementors do not have and made every client null-check a raw pointer.
//  The polarized wave function is the PRIMARY type here, not a special case bolted onto the base --
//  an unpolarized one simply does not answer this question, and now cannot be asked it.
//
//  The idiom is the one qcChargeDensity already uses for exactly this shape (tSpinResolved_CD): a
//  data-free face that is NOT a tWaveFunction, reached by the sanctioned abstract->abstract
//  dynamic_cast, so capabilities live only on the types that have them.
//
export template <class T> class tSpinResolvedWF
{
public:
    typedef ScalarFunction<double> sf_t;
    virtual ~tSpinResolvedWF() {}
    //! \brief BUILDS \f$m(r)\f$ over BOTH channels -- ALLOCATES, hence the owning return (V1.25).
    virtual std::unique_ptr<sf_t> GetSpinDensity() const=0;
};

export using SpinResolvedWF  = tSpinResolvedWF<double>;
export using cSpinResolvedWF = tSpinResolvedWF<dcmplx>;

} //namespace

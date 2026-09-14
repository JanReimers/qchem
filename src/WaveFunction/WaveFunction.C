// File: WaveFunction.C  Interface for a wave function.
module;
#include <memory>   // std::unique_ptr -- GetChargeDensity BUILDS its density (V1.25)
#include <vector>
export module qchem.WaveFunction;
export import qchem.EnergyLevel;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;   // Irrep + Spin + SpinGroup (the imposed spin subgroup)
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
    //! \brief THE IMPOSED SPIN SUBGROUP this state was built under (V1.37) -- a property of the run, like
    //! the point group, not of a type: \c UnPolarized = one folded doublet per spatial irrep (Spin::None),
    //! \c Polarized = the collinear Up/Down pair.
    virtual SpinGroup       GetSpinGroup    () const=0;
    //! \brief BUILDS the collinear magnetization \f$m(r)=\rho_\uparrow-\rho_\downarrow\f$ -- ALLOCATES,
    //! hence the owning return (V1.25).  Under imposed SU(2) (\c UnPolarized) \f$m\equiv0\f$ by symmetry
    //! and this THROWS rather than raster a field of zeros: the caller has \c GetSpinGroup() to ask first,
    //! which is exactly what the facade does.  (This used to be a separate \c tSpinResolvedWF capability
    //! face carried by the polarized TYPE; with one wave-function class the subgroup is the only thing
    //! left to ask, so it is asked.)
    virtual std::unique_ptr<sf_t> GetSpinDensity() const=0;
    // (No Emit*() here.  V1.14: a class does not tell another class WHEN to report -- the composite WF
    // announces its basis usage itself, from FillOrbitals, the moment the occupations exist.)

private:
    tWaveFunction& operator=(const tWaveFunction&);
};

export using WaveFunction  = tWaveFunction<double>;
export using cWaveFunction = tWaveFunction<dcmplx>;


} //namespace

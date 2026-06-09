// File: WaveFunction.C  Interface for a wave function.
module;
#include <vector>
export module qchem.WaveFunction;
export import qchem.EnergyLevel;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;
export import qchem.Symmetry.ElectronConfiguration;
import qchem.ScalarFunction;
export import qchem.Orbitals;


namespace qchem::WaveFunction
{

export using qchem::ChargeDensity::DM_CD;
export using qchem::Orbitals::EnergyLevels;
using Hamiltonian::Hamiltonian;
using Orbitals::Orbitals; //Keep this one last, otherwise it interferes with the two previous declarations!

export class WaveFunction 
{
public:
    typedef ScalarFunction<double> sf_t;
    typedef std::vector<Irrep> iqns_t;
    virtual ~WaveFunction() {};

    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*)      =0;
    virtual const Orbitals* GetOrbitals     (const Irrep&         ) const=0;
    virtual       Orbitals* GetOrbitals     (const Irrep&         )      =0;
    virtual void            FillOrbitals    (double mergeTol)=0; //WF knows internally the electronic structure
    virtual DM_CD*          GetChargeDensity() const=0;
    virtual sf_t*           GetSpinDensity  () const=0; //Returns a null ptr for un polarized WF.
    virtual EnergyLevels    GetEnergyLevels () const=0; 
    virtual iqns_t          GetQNs          () const=0;
    virtual void            DisplayEigen    () const=0;


private:
    WaveFunction& operator=(const WaveFunction&);
};

} //namespace

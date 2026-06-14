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
    // Direct-minimization hooks (cf. SCFIterator direct-min loop):
    //   build the Fock and ask each accelerator to compute its step (no orbital move);
    //   returns false in the seed step (the caller should DoSCFIteration to diagonalize).
    virtual bool            BuildFockAndComputeSteps(Hamiltonian&,const DM_CD*) {return false;}
    //   move the orbitals to geodesic fraction t (commit=false is a line-search trial) and refill.
    virtual void            MoveOrbitals    (double t, bool commit, double mergeTol) {}
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

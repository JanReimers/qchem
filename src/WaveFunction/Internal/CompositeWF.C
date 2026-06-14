// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <vector>
#include <map>
#include <memory>
export module qchem.WaveFunction.Internal.CompositeWF;
export import qchem.WaveFunction.SCF;
import qchem.SCFAccelerator;
import qchem.WaveFunction.Internal.IrrepWF;
import qchem.WaveFunction.Types;

export namespace qchem::WaveFunction
{

using SCFAccelerators::SCFAccelerator;

class CompositeWF
    : public virtual SCFWaveFunction
{
public:
    CompositeWF(const bs_t*,const ElectronConfiguration*,SCFAccelerator*);
    ~CompositeWF();

    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual bool            BuildFockAndComputeSteps(Hamiltonian&,const DM_CD*);
    virtual void            MoveOrbitals    (double t, bool commit, double mergeTol);
    virtual const Orbitals* GetOrbitals     (const Irrep&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;} 
    virtual void            FillOrbitals    (double mergeTol);
    virtual iqns_t          GetQNs          () const;

    virtual DM_CD*          GetChargeDensity(Spin) const;
    virtual EnergyLevels    GetEnergyLevels (Spin) const; 

    
protected:
    void MakeIrrepWFs(Spin);

private:
    typedef std::unique_ptr<IrrepWF> uiwf_t;

    const bs_t*                  itsBS; 
    const ElectronConfiguration* itsEC;
    SCFAccelerator*              itsAccelerator;
    EnergyLevels                 itsELevels;
    LAParams                     itsLAParams; //Numerical control of general eigen solution.

    std::map<Spin,EnergyLevels>  itsSpin_ELevels;

    std::vector<uiwf_t>                   itsIWFs;
    std::map<Irrep,IrrepWF*>         itsQNWFs; //sort by Irrep for easy lookup.
    std::map<Spin,std::vector<IrrepWF*>> itsSpinWFs; //Sort by spin.
};

} //namespace

// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <vector>
#include <map>
#include <memory>
import qchem.SCFAccelerator;
export module qchem.WaveFunction.Internal.CompositeWF;
export import qchem.WaveFunction;
export import qchem.BasisSet;
import qchem.WaveFunction.Internal.IrrepWF;

export class CompositeWF
    : public virtual WaveFunction
{
public:
    CompositeWF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator*);
    ~CompositeWF();

    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual const Orbitals* GetOrbitals     (const Irrep_QNs&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep_QNs&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;} 
    virtual void            FillOrbitals    (double mergeTol);
    virtual iqns_t          GetQNs          () const;

    virtual DM_CD*          GetChargeDensity(Spin) const;
    virtual EnergyLevels    GetEnergyLevels (Spin) const; 

protected:
    void MakeIrrepWFs(Spin);

private:
    typedef std::unique_ptr<IrrepWF> uiwf_t;

    const BasisSet*              itsBS; 
    const ElectronConfiguration* itsEC;
    SCFAccelerator*              itsAccelerator;
    EnergyLevels                 itsELevels;
    std::map<Spin,EnergyLevels>  itsSpin_ELevels;

    std::vector<uiwf_t>                   itsIWFs;
    std::map<Irrep_QNs,IrrepWF*>         itsQNWFs; //sort by Irrep for easy lookup.
    std::map<Spin,std::vector<IrrepWF*>> itsSpinWFs; //Sort by spin.
};


// File: Composite_WF.H  Wave function as a list of Irrep wave functions.
#ifndef _Composite_WF_H_
#define _Composite_WF_H_

#include <WaveFunction/WaveFunction.H>
#include <EnergyLevel.H>
#include <vector>
#include <map>
#include <memory>

class Irrep_WF;

class Composite_WF
    : public virtual WaveFunction
{
public:
    Composite_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator*);
    ~Composite_WF();

    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual const Orbitals* GetOrbitals     (const Irrep_QNs&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep_QNs&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;} 
    virtual void            FillOrbitals    (double mergeTol);
    virtual iqns_t          GetQNs          () const;

    virtual DM_CD*          GetChargeDensity(Spin) const;
    virtual EnergyLevels    GetEnergyLevels (Spin) const; 

protected:
    void MakeIrrep_WFs(Spin);

private:
    typedef std::unique_ptr<Irrep_WF> uiwf_t;

    const BasisSet*              itsBS; 
    const ElectronConfiguration* itsEC;
    SCFAccelerator*              itsAccelerator;
    EnergyLevels                 itsELevels;
    std::map<Spin,EnergyLevels>  itsSpin_ELevels;

    std::vector<uiwf_t>                   itsIWFs;
    std::map<Irrep_QNs,Irrep_WF*>         itsQN_WFs; //sort by Irrep for easy lookup.
    std::map<Spin,std::vector<Irrep_WF*>> itsSpin_WFs; //Sort by spin.
};

#endif //_Composite_WF_H_

// File: Polarized_WF.H  Wave function for an polarized system.
#ifndef _Polarized_WF_H_
#define _Polarized_WF_H_


#include "WaveFunction.H"
#include <EnergyLevel.H>
#include <vector>
#include <map>
#include <memory>

class Polarized_WF
    : public virtual WaveFunction
{
public:
    Polarized_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator& acc);
    
    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual const Orbitals* GetOrbitals     (const Irrep_QNs&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep_QNs&)      ;
    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const; 
    virtual EnergyLevels    GetEnergyLevels () const; 
    virtual void            FillOrbitals    ();
    virtual iqns_t          GetQNs          () const;
    virtual void            DisplayEigen    () const;

private:
    typedef std::unique_ptr<Irrep_WF> uiwf_t;
    const BasisSet* itsBS; 
    const ElectronConfiguration* itsEC;
    EnergyLevels itsUpELevels;
    EnergyLevels itsDnELevels;
    std::vector<uiwf_t> itsSpinUpIWFs;
    std::vector<uiwf_t> itsSpinDnIWFs;
    std::map<Irrep_QNs,Irrep_WF*> itsQN_WFs;
};

#endif //_Polarized_WF_H_

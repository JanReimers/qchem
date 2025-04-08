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
    Polarized_WF(const BasisSet*,const ElectronConfiguration*);
    
    virtual void             DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual DM_CD*           GetChargeDensity() const;
    virtual const Orbitals*  GetOrbitals     (const Irrep_QNs&) const;
    virtual void             FillOrbitals    (const ElectronConfiguration*);
    virtual void             DisplayEigen    (               ) const;

private:
    typedef std::unique_ptr<Irrep_WF> uiwf_t;
    const BasisSet* itsBS; 
    const ElectronConfiguration* itsEC;
    EnergyLevels itsUpELevels;
    EnergyLevels itsDnELevels;
    std::vector<uiwf_t> itsSpinUpIWFs;
    std::vector<uiwf_t> itsSpinDnIWFs;
    std::map<Irrep_QNs,const Irrep_WF*> itsQN_WFs;
};

#endif //_Polarized_WF_H_

// File: UnPolarized_WF.H  Wave function for an unpolarized atom.
#ifndef _UnPolarized_WF_H_
#define _UnPolarized_WF_H_



#include "WaveFunction.H"
#include <EnergyLevel.H>
#include <vector>
#include <map>
#include <memory>

class Irrep_WF;

class UnPolarized_WF
    : public virtual WaveFunction
{
public:
    UnPolarized_WF();
    UnPolarized_WF(const BasisSet*,const ElectronConfiguration*,SCFAccelerator& acc);
    virtual ~UnPolarized_WF();

    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual const Orbitals* GetOrbitals     (const Irrep_QNs&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep_QNs&)      ;
    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;} 
    virtual void            FillOrbitals    ();
    virtual iqns_t          GetQNs          () const;
    virtual void            DisplayEigen    () const;

private:
    typedef std::unique_ptr<Irrep_WF> uiwf_t;
    const BasisSet* itsBS; 
    const ElectronConfiguration* itsEC;
    EnergyLevels itsELevels;
    std::vector<uiwf_t> itsIWFs;
    std::map<Irrep_QNs,Irrep_WF*> itsQN_WFs;
};

#endif //_UnPolarized_WF_H_

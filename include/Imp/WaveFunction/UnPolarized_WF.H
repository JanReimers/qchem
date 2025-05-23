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
    UnPolarized_WF(const BasisSet*,const ElectronConfiguration*);
    virtual ~UnPolarized_WF();

    virtual void            DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    virtual DM_CD*          GetChargeDensity() const;
    virtual sf_t*           GetSpinDensity  () const {return 0;}
    virtual const Orbitals* GetOrbitals     (const Irrep_QNs&) const;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;} 
    virtual void            FillOrbitals    ();
    virtual void            DisplayEigen    (               ) const;

private:
    typedef std::unique_ptr<Irrep_WF> uiwf_t;
    const BasisSet* itsBS; 
    const ElectronConfiguration* itsEC;
    EnergyLevels itsELevels;
    std::vector<uiwf_t> itsIWFs;
    std::map<Irrep_QNs,const Irrep_WF*> itsQN_WFs;
};

#endif //_UnPolarized_WF_H_

// File: MasterPolarizedWF.H  Wave function for an polarized system.
#ifndef _MasterPolarizedWF_H_
#define _MasterPolarizedWF_H_


#include "WaveFunction.H"
#include <EnergyLevel.H>
#include <vector>
#include <map>
#include <memory>

class IrrepWaveFunction;

class MasterPolarizedWF
    : public virtual WaveFunction
{
public:
    MasterPolarizedWF(const BasisSet*,const ElectronConfiguration*);
    
    virtual void          DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*     GetChargeDensity() const;
    virtual Orbitals*     GetOrbitals     (const Irrep_QNs&) const;
    virtual void          FillOrbitals    (const ElectronConfiguration*);
    virtual void          DisplayEigen    (               ) const;

private:
    typedef std::unique_ptr<IrrepWaveFunction> uiwf_t;
    const BasisSet* itsBS; 
    const ElectronConfiguration* itsEC;
    EnergyLevels itsUpELevels;
    EnergyLevels itsDnELevels;
    std::vector<uiwf_t> itsSpinUpIWFs;
    std::vector<uiwf_t> itsSpinDnIWFs;
    std::map<Irrep_QNs,const IrrepWaveFunction*> itsQN_WFs;
};

#endif //_MasterPolarizedWF_H_

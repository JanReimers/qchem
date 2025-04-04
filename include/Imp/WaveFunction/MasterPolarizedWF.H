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
    : public virtual PolarizedWF
{
public:
    MasterPolarizedWF(const BasisSet*,const ElectronConfiguration*);
    
    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*           GetChargeDensity(Spin           ) const;
    virtual Orbitals*           GetOrbitals     (const Irrep_QNs&) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double nElectrons);
    virtual void                DisplayEigen    (               ) const;

    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

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

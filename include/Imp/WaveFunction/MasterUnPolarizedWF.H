// File: MasterUnPolarizedWF.H  Wave function for an unpolarized atom.
#ifndef _MasterUnPolarizedWF_H_
#define _MasterUnPolarizedWF_H_



#include "WaveFunction.H"
#include <EnergyLevel.H>
#include <vector>
#include <map>
#include <memory>

class IrrepWaveFunction;

class MasterUnPolarizedWF
    : public virtual WaveFunction
{
public:
    MasterUnPolarizedWF();
    MasterUnPolarizedWF(const BasisSet*,const ElectronConfiguration*);
    virtual ~MasterUnPolarizedWF();

    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*           GetChargeDensity(  ) const;
    virtual Orbitals*           GetOrbitals     (const Irrep_QNs&) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double nElectrons);
    virtual void                DisplayEigen    (               ) const;

private:
    typedef std::unique_ptr<IrrepWaveFunction> uiwf_t;
    const BasisSet* itsBS; 
    const ElectronConfiguration* itsEC;
    EnergyLevels itsELevels;
    std::vector<uiwf_t> itsIWFs;
    std::map<Irrep_QNs,const IrrepWaveFunction*> itsQN_WFs;
};

#endif //_MasterUnPolarizedWF_H_

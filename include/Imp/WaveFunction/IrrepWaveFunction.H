// File: IrrepWaveFunction.H  Wave function for an irreducable representation.
#ifndef _IrrepWaveFunction_H_
#define _IrrepWaveFunction_H_



#include <WaveFunction.H>
#include <Spin.H>
#include <EnergyLevel.H>
#include <Orbital_QNs.H>

class IrrepWaveFunction
: public virtual PMStreamableObject
{
public:
    IrrepWaveFunction(const TOrbital_IBS<double>*, const Spin&);
    virtual ~IrrepWaveFunction();

    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*           GetChargeDensity(Spin           ) const;
    virtual Orbitals*           GetOrbitals     (const Irrep_QNs&) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    // virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double nElectrons);
    virtual void                DisplayEigen    (               ) const;

    virtual std::ostream&      Write(std::ostream&) const;
    virtual std::istream&      Read (std::istream&)      ;
    static  IrrepWaveFunction* Factory(std::istream&)        ;
private:
    IrrepWaveFunction(const IrrepWaveFunction&);

    Orbitals*     itsOrbitals;
    Spin          itsSpin;
    const Symmetry*    itsQN;
    Irrep_QNs     itsQNs;
    EnergyLevels  itsELevels;

}; 

#endif //IrrepWaveFunction

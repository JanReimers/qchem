// File: IrrepWaveFunction.H  Wave function for an irreducable representation.
#ifndef _IrrepWaveFunction_H_
#define _IrrepWaveFunction_H_



#include <WaveFunction.H>
#include <Spin.H>
#include <EnergyLevel.H>
#include <Orbital_QNs.H>

class IrrepWaveFunction
    : public virtual WaveFunction
{
public:
    IrrepWaveFunction(const TOrbital_IBS<double>*, const Spin&);
    virtual ~IrrepWaveFunction();

    virtual void                DoSCFIteration  (Hamiltonian&   )      ;
    virtual Exact_CD*           GetChargeDensity(Spin           ) const;
    virtual Orbitals*           GetOrbitals     (const QNs&,Spin) const;
    virtual const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    virtual SCFIterator*        MakeIterator    (Hamiltonian*, Exact_CD* guess, double nElectrons);
    virtual void                DisplayEigen    (               ) const;

    virtual std::ostream&      Write(std::ostream&) const;
    virtual std::istream&      Read (std::istream&)      ;

private:
    IrrepWaveFunction(const IrrepWaveFunction&);

    Orbitals*     itsOrbitals;
    Spin          itsSpin;
    const QNs*    itsQN;
    Orbital_QNs   itsQNs;
    EnergyLevels  itsELevels;

}; 

#endif //IrrepWaveFunction

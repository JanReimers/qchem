// File: IrrepWaveFunction.H  Wave function for an irreducable representation.
#ifndef _IrrepWaveFunction_H_
#define _IrrepWaveFunction_H_



#include "WaveFunction.H"
#include "BasisSet.H"
#include "Misc/Spin.H"

class IrrepWaveFunction
    : public virtual WaveFunction
{
public:
    IrrepWaveFunction();
    IrrepWaveFunction(const IrrepBasisSet*, const Spin&);
    virtual ~IrrepWaveFunction();

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           UpdateElectronDumper(ElectronDumper&)      ;
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double nElectrons, double kT, bool showplot);

    virtual std::ostream&      Write(std::ostream&) const;
    virtual std::istream&      Read (std::istream&)      ;

private:
    IrrepWaveFunction(const IrrepWaveFunction&);

    Orbitals*            itsOrbitals;
    Spin                 itsSpin;
}; 

#endif //IrrepWaveFunction

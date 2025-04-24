// File: Irrep_WF.H  Wave function for an irreducable representation.
#ifndef _IrrepWaveFunction_H_
#define _IrrepWaveFunction_H_

#include <EnergyLevel.H>
#include <Orbital_QNs.H>

class Irrep_WF
{
public:
    Irrep_WF(const TOrbital_IBS<double>*, const Spin&);
    ~Irrep_WF();

    void                DoSCFIteration  (Hamiltonian&,const DM_CD*   )      ;
    DM_CD*              GetChargeDensity(                            ) const;
    const Orbitals*     GetOrbitals     (           ) const;
    const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    void                DisplayEigen    (                            ) const;
    const Irrep_QNs&    GetQNs          () const {return itsIrrep;}
 private:
    Irrep_WF(const Irrep_WF&);

    Orbitals*     itsOrbitals; //Owned
    Irrep_QNs     itsIrrep;
    EnergyLevels  itsELevels;

}; 

#endif //IrrepWaveFunction

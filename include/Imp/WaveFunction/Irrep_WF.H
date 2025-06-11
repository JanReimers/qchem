// File: Irrep_WF.H  Wave function for an irreducable representation.
#ifndef _IrrepWaveFunction_H_
#define _IrrepWaveFunction_H_

#include <EnergyLevel.H>
#include <Orbital_QNs.H>
#include "oml/smatrix.h"

class Irrep_WF
{
public:
    Irrep_WF(const TOrbital_IBS<double>*, LASolver<double>*, const Irrep_QNs& ,SCFIrrepAccelerator*);
    ~Irrep_WF();

    void                CalculateH      (Hamiltonian&,const DM_CD*   )      ;
    void                DoSCFIteration  ()      ;
    DM_CD*              GetChargeDensity() const;
    const Orbitals*     GetOrbitals     () const;
          Orbitals*     GetOrbitals     ()      ;
    const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    void                DisplayEigen    () const;
    const Irrep_QNs&    GetQNs          () const {return itsIrrep;}
    Vector<double>      Get_BS_Diagonal () const;

 private:
    Irrep_WF(const Irrep_WF&);

    const TOrbital_IBS<double>*  itsBasisSet;
    LASolver<double>*            itsLASolver;
    TOrbitals<double>*           itsOrbitals; //Owned
    Irrep_QNs                    itsIrrep;
    EnergyLevels                 itsELevels;
    SCFIrrepAccelerator*         itsAccelerator;
    SMatrix<double>              itsDPrime,itsF; // DPrime=C'*Cd',  U*D*Ud, D=C*Cd (outer product)
}; 

#endif //IrrepWaveFunction

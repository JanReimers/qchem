// File: IrrepWF.H  Wave function for an irreducable representation.
module;
export module qchem.WaveFunction.Internal.IrrepWF;
import qchem.SCFAccelerator;
export import qchem.EnergyLevel;
export import qchem.LASolver;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.IrrepBasisSet;
export import qchem.Symmetry.Irrep;
export import qchem.Symmetry.ElectronConfiguration;

import oml;

export class IrrepWF
{
public:
    IrrepWF(const Real_OIBS*, LASolver<double>*, const Irrep_QNs& ,SCFIrrepAccelerator*);
    ~IrrepWF();

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
    IrrepWF(const IrrepWF&);

    const Real_OIBS*  itsBasisSet;
    LASolver<double>*            itsLASolver;
    TOrbitals<double>*           itsOrbitals; //Owned
    Irrep_QNs                    itsIrrep;
    EnergyLevels                 itsELevels;
    SCFIrrepAccelerator*         itsAccelerator;
    SMatrix<double>              itsDPrime,itsF; // DPrime=C'*Cd',  U*D*Ud, D=C*Cd (outer product)
}; 


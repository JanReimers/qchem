// File: IrrepWF.H  Wave function for an irreducable representation.
module;
export module qchem.WaveFunction.Internal.IrrepWF;
import qchem.WaveFunction.Types;
import qchem.SCFAccelerator;
export import qchem.EnergyLevel;
export import qchem.Orbitals;
export import qchem.LASolver;
export import qchem.Hamiltonian;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Irrep;
export import qchem.Symmetry.ElectronConfiguration;

export namespace qchem::WaveFunction
{

using ChargeDensity::DM_CD;
using Hamiltonian::Hamiltonian;
using Orbitals::TOrbitals;
using Orbitals::EnergyLevels;
using Orbitals::Orbitals; //Keep this one here, otherwise it interferes with the two previous declarations!
using SCFAccelerators::SCFIrrepAccelerator;

class IrrepWF
{
public:
    IrrepWF(const obs_t*, LASolver<double>*, const Irrep& ,SCFIrrepAccelerator*);
    ~IrrepWF();

    void                CalculateH      (Hamiltonian&,const DM_CD*   )      ;
    void                DoSCFIteration  ()      ;
    bool                ComputeStep     ()      ; //direct-min: accelerator computes its step
    void                MoveOrbitals    (double t, bool commit)      ; //move to geodesic fraction t
    DM_CD*              GetChargeDensity() const;
    const Orbitals*     GetOrbitals     () const;
          Orbitals*     GetOrbitals     ()      ;
    const EnergyLevels& FillOrbitals    (const ElectronConfiguration*);
    void                DisplayEigen    () const;
    const Irrep&    GetQNs          () const {return itsIrrep;}
    rvec_t      Get_BS_Diagonal () const;

 private:
    IrrepWF(const IrrepWF&);

    const obs_t*         itsBasisSet;
    LASolver<double>*    itsLASolver;
    TOrbitals<double>*   itsOrbitals; //Owned
    Irrep            itsIrrep;
    EnergyLevels         itsELevels;
    SCFIrrepAccelerator* itsAccelerator;
    rsmat_t              itsDPrime; // DPrime=C'*Cd',  U*D*Ud, D=C*Cd (outer product)
    rsmat_t              itsF;
}; 

} //namespace

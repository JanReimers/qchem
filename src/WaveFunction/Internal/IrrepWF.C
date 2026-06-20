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
export import qchem.ElectronConfiguration;

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
    const EnergyLevels& FillOrbitals    (double ne); //occupy with a given electron count (aufbau)

    // Maximum Overlap Method (MOM): score each *current* orbital by how much it overlaps the
    // reference occupied subspace (previous iteration's occupied orbitals), so occupation can
    // follow orbital character instead of eigenvalue.  Empty if no reference captured yet.
    rvec_t              MOMScores       () const;
    void                CaptureMOMReference()      ; //snapshot the occupied orbitals as the next reference

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
    rmat_t               itsRefOccCPrime; // MOM reference: occupied C' columns (nbasis x nocc); empty=none
};

} //namespace

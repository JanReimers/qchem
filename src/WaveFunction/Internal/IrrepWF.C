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

using ChargeDensity::rDM_CD;
using ChargeDensity::tDM_CD;
using ChargeDensity::tChargeDensity;
using qchem::Hamiltonian::tHamiltonian;
// tbs_t (the whole composite basis = the cross-irrep view threaded to the dynamic Fock terms) comes from
// qchem.WaveFunction.Types; it is the same BasisSet::tBasisSet<T> the Hamiltonian layer's GetMatrix expects.
using Orbitals::TOrbitals;
using Orbitals::EnergyLevels;
using Orbitals::Orbitals; //Keep this one here, otherwise it interferes with the two previous declarations!
using SCFAccelerators::tSCFIrrepAccelerator;

// Per-irrep wave function, templated on the matrix element type T (rX/cX); IrrepWF is the <double>
// alias (atom/molecule), cIrrepWF the dcmplx (plane-wave / Bloch-irrep) instantiation.
template <class T> class tIrrepWF
{
public:
    tIrrepWF(const tobs_t<T>*, LASolver<T>*, const Irrep& ,tSCFIrrepAccelerator<T>*);
    ~tIrrepWF();

    void                CalculateH      (tHamiltonian<T>&,const tChargeDensity<T>*,const tbs_t<T>* wholeBasis);
    void                DoSCFIteration  ()      ;
    bool                ComputeStep     ()      ; //direct-min: accelerator computes its step
    void                MoveOrbitals    (double t, bool commit)      ; //move to geodesic fraction t
    tDM_CD<T>*          GetChargeDensity() const;
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
    const Irrep&    GetIrrep        () const {return itsIrrep;}   // this WF's irrep (the proper map key)
    rvec_t      Get_BS_Diagonal () const;

 private:
    tIrrepWF(const tIrrepWF&);

    const tobs_t<T>*         itsBasisSet;
    LASolver<T>*             itsLASolver;
    TOrbitals<T>*            itsOrbitals; //Owned
    Irrep                    itsIrrep;
    EnergyLevels             itsELevels;
    tSCFIrrepAccelerator<T>* itsAccelerator;
    hmat_t<T>                itsDPrime; // DPrime=C'*Cd',  U*D*Ud, D=C*Cd (outer product)
    hmat_t<T>                itsF;
    mat_t<T>                 itsRefOccCPrime; // MOM reference: occupied C' columns (nbasis x nocc); empty=none
    int                      itsFillCount=0;  // # of FillOrbitals calls (≈ SCF iteration) -- IMOM capture delay
};

using IrrepWF  = tIrrepWF<double>;
using cIrrepWF = tIrrepWF<dcmplx>;

} //namespace

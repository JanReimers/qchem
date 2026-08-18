// File: IrrepWF.H  Wave function for an irreducable representation.
module;
#include <memory>   // std::unique_ptr -- GetChargeDensity BUILDS its density (V1.25)
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
export import qchem.ElectronConfiguration.OccupationPolicy;  // the fills consult the iterator's policy (V1.11 inc 3)

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
    //! Step 3c-2: the REAL child's Fock build inside a COMPLEX run -- drives the Hamiltonian's
    //! Ham_RealBlock assembly face with the run's (complex) density and whole-basis view; itsF comes back
    //! real.  Declared for every T (the <dcmplx> instantiation's body throws -- a complex child never
    //! consumes this face); MEANINGFUL only on tIrrepWF<double>.
    void                CalculateH      (qchem::Hamiltonian::Ham_RealBlock&,
                                         const tChargeDensity<dcmplx>*,const tbs_t<dcmplx>* wholeBasis);
    void                DoSCFIteration  ()      ;
    bool                ComputeStep     ()      ; //direct-min: accelerator computes its step
    void                MoveOrbitals    (double t, bool commit)      ; //move to geodesic fraction t
    std::unique_ptr<tDM_CD<T>> GetChargeDensity() const;   //!< BUILDS it (V1.25: owning return)
    const Orbitals*     GetOrbitals     () const;
          Orbitals*     GetOrbitals     ()      ;
    //! Occupy with a given electron count: the POLICY decides the fill (DecideBlockFill -- occupancy rule +
    //! ranking, including the direct minimiser's HeldOccupationPolicy keeping the stored block), the
    //! orbitals execute, and this block accumulates its BZ-weighted −TS into \a pol (V1.11).
    const EnergyLevels& FillOrbitals    (OccupationPolicy<T>& pol, double ne);
    //! Occupy at a GIVEN chemical potential μ (the global-μ metal fill, doc/GPWPlan1.md item 3): empties, sets
    //! g_i·f_i on THIS block's orbitals at \a mu (plain energy Fermi -- the composite solved μ on bare ε across
    //! the mesh), stores this block's D'/levels and accumulates its −TS into \a pol.  Requires smearing on.
    //! No MOM (a metal fills by energy).
    const EnergyLevels& FillOrbitalsAtMu(OccupationPolicy<T>& pol, double mu);
    // NB: the MOM machinery (scores, reference capture/adopt/release) and the smearing/MOM configuration
    // are GONE from this class -- they are OccupationPolicy state, keyed by this block's Irrep (V1.11
    // inc 3).  The −TS side-channel (GetEntropyTerm) went with them: each fill accumulates
    // w_k·(−TS_block) into the policy, and the SCFIterator reads the run total from there.

    void                DisplayEigen    () const;
    const Irrep&    GetIrrep        () const {return itsIrrep;}   // this WF's irrep (the proper map key)
    rvec_t      Get_BS_Diagonal () const;
    //! Occupation-weighted Mulliken gross population per basis function (the basis-usage heat map, doc/GPWPlan1
    //! §1): delegate to the orbitals (which own D) with THIS irrep's AO overlap.  Sum = the irrep's electrons.
    rvec_t      GetBasisPopulations() const {return itsOrbitals->GetBasisPopulations(itsBasisSet->Overlap());}

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
    // (The MOM reference / fill-count / smearing-config members moved to OccupationPolicy -- V1.11 inc 3.)
};

using IrrepWF  = tIrrepWF<double>;
using cIrrepWF = tIrrepWF<dcmplx>;

} //namespace

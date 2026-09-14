// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <vector>
#include <map>
#include <memory>
#include <variant>
export module qchem.WaveFunction.Internal.CompositeWF;
export import qchem.WaveFunction.SCF;
import qchem.SCFAccelerator;
import qchem.WaveFunction.Internal.IrrepWF;
import qchem.WaveFunction.Types;
import qchem.LASolver;   // qchem::Ortho (the basis-overlap orthogonalisation: Cholesky | Eigen | SVD + tol)
import qchem.Matrix3D;   // Matrix3D -- the reciprocal point ops passed to each composite density (IBZ symmetrization)

export namespace qchem::WaveFunction
{

using SCFAccelerators::SCFAccelerator;
using qchem::Hamiltonian::tHamiltonian;
using ChargeDensity::tDM_CD;
using ChargeDensity::tChargeDensity;

//! \brief THE WF CHILD SLOT (doc/RealComplexPlan.md §4, Step 2b): one per-irrep wave function, typed by
//! ITS block's scalar rather than the composite's face -- the WF twin of tComposite_CD's cd_child_t, so a
//! real TRIM block's eigensolve/C/D can live beside general-k complex blocks on a mixed mesh.  Aggregation
//! stays single-source: scalar-signature calls (GetIrrep, GetOrbitals -- the Orbitals base is non-template,
//! DoSCFIteration, ComputeStep, MoveOrbitals, GetChargeDensity -- 2a's typed Inserts take either) visit with
//! ONE generic lambda; the T-typed calls (CalculateH's Hamiltonian, the fills' OccupationPolicy) forward to
//! the same-scalar alternative, with the cross arm throwing until Step 2c/3 land the narrowing.
using iwf_child_t = std::variant<std::unique_ptr<tIrrepWF<double>>, std::unique_ptr<tIrrepWF<dcmplx>>>;
//! Non-owning mirror of the child slot (the Irrep/Spin lookup maps + the reservoir grouping).
using iwf_ref_t   = std::variant<tIrrepWF<double>*, tIrrepWF<dcmplx>*>;

// THE wave function: ONE composite of per-irrep wave functions over the FULL Irrep (spatial ⊗ Spin) --
// doc/CleanupCandidates.md V1.37.  Templated on the matrix element type T (rX/cX); CompositeWF is the
// <double> alias (atoms/molecules), cCompositeWF the <dcmplx> (plane-wave / single-k Bloch-irrep)
// instantiation.
//
// Pol/UnPol is the IMPOSED SPIN SUBGROUP (qchem::SpinGroup), a ctor argument like the point group is a
// property of the basis -- NOT a pair of thin subclasses (tPolarizedWF / tUnPolarizedWF are gone).  It
// decides only which spin irreps the per-irrep children are built for: Spin::None (one folded doublet
// per spatial irrep, degeneracy 2 -- the CLAUDE.md "UnPol is the efficient special case", carried by the
// LABEL) or Spin::Up + Spin::Down.  Everything else -- the fills, the density, the levels -- is one code
// path over the children.  The subgroup shows in exactly two places: the level DISPLAY (a side-by-side
// ↑/↓ table vs. one column) and GetSpinDensity (m ≡ 0 under SU(2), so it is not built).
template <class T> class tCompositeWF
    : public virtual tSCFWaveFunction<T>
{
public:
    typedef typename tWaveFunction<T>::iqns_t iqns_t;
    typedef typename tWaveFunction<T>::sf_t   sf_t;

    //! \a g is the imposed spin subgroup (see the class note).  \a basisOrtho selects how the (per-irrep)
    //! orbital-overlap S is orthogonalised for the generalised eigenproblem: \c Cholesky (default; requires S
    //! positive-definite) or \c Eigen / \c SVD with a \a basisOrthoTol cutoff that DROPS near-null
    //! eigen/singular values -- canonical orthogonalisation for a linearly-dependent basis (e.g. diffuse
    //! Gaussians on a dense lattice).  \a basisOrthoTol\f$\le0\f$ = keep all.
    tCompositeWF(const tbs_t<T>*,const ElectronConfiguration*,SpinGroup g,SCFAccelerator*,
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0);
    ~tCompositeWF();

    virtual void            DoSCFIteration  (tHamiltonian<T>&,const tChargeDensity<T>*   )      ;
    //! Iteration-0 seed.  BUILDS and hands over the first real density (V1.25: owning return).
    virtual std::unique_ptr<tDM_CD<T>> Init(tHamiltonian<T>&,const tChargeDensity<T>*,
                                            OccupationPolicy<T>&, double mergeTol);
    virtual bool            BuildFockAndComputeSteps(tHamiltonian<T>&,const tChargeDensity<T>*);
    virtual void            MoveOrbitals    (OccupationPolicy<T>&, double t, bool commit, double mergeTol);
    virtual const Orbitals* GetOrbitals     (const Irrep&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;}   //!< all spin irreps merged
    virtual void            FillOrbitals    (OccupationPolicy<T>&, double mergeTol);
    // (SetMOM/SetSmearing/GetEntropyTerm/AdoptMOMReference/ReleaseMOMReference are GONE -- the
    //  SCFIterator's OccupationPolicy slot owns that configuration and state, V1.11 inc 3.)
    virtual iqns_t          GetQNs          () const;
    virtual void            DisplayEigen    () const;
    virtual SpinGroup       GetSpinGroup    () const {return itsSpinGroup;}

    //! BUILDS the whole-system density (V1.25): ONE tComposite_CD over EVERY child block, in the order
    //! they were built (spin irrep by spin irrep) -- a polarized run's Up blocks then its Down blocks, and
    //! the composite's channel views are what a spin-native consumer reads.
    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity() const;
    //! BUILDS one spin irrep's density alone (V1.25) -- the channel as a composite of its own, for the
    //! spin density and for tests; the SCF consumes the whole-system one above.
    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity(Spin) const;
    //! \copydoc tWaveFunction::GetSpinDensity
    virtual std::unique_ptr<sf_t> GetSpinDensity() const;
    virtual EnergyLevels    GetEnergyLevels (Spin) const;


protected:
    void MakeIrrepWFs(Spin);
    //! One block's WF child, typed by the BLOCK's scalar U (Step 3c-2): the single-source body behind the
    //! mixed MakeIrrepWFs walk -- U==T is the native path, U==double under a complex face is the real
    //! TRIM child (its accelerator comes from the §6 typed Create; its density Inserts 2a's <double> arm).
    template <class U> void MakeOneIrrepWF(const tobs_t<U>*, Spin);

private:
    //! Announce the run report's `basis.usage` block (per-function occupation-weighted populations) -- called
    //! by FillOrbitals ITSELF, at the moment the occupations exist (V1.14: providers self-report at their
    //! own trigger; nobody tells this class when).  EmitAt is idempotent and run-scoped, so announcing on
    //! every fill is cheap and the json ends holding the LAST (converged) fill.  No-op when no run is open.
    void AnnounceBasisUsage() const;
    typedef tIrrepWF<T> iwf_t;   // the same-face child kind (what MakeIrrepWFs builds today)
    //! RANKED integer fill of one reservoir (the molecular cross-irrep aufbau, one spin channel): pick which
    //! orbitals across the reservoir's blocks are occupied, then fill each block with its resulting count.
    void FillReservoirRanked    (OccupationPolicy<T>&, const std::vector<iwf_ref_t>&, double mergeTol, bool useMOM);
    //! Smeared fill of one reservoir: solve ONE μ over the reservoir's blocks -- across the Bloch mesh when
    //! the partition spans spatial (doc/GPWPlan1.md item 3), across SPIN when it spans spin (free moment).
    void FillReservoirAtSharedMu(OccupationPolicy<T>&, const std::vector<iwf_ref_t>&, double mergeTol);
    //! The two level tables: the subgroup picks one (see the class note).
    void DisplayEigenPolarized  () const;
    void DisplayEigenUnPolarized() const;

    const tbs_t<T>*              itsBS;
    const ElectronConfiguration* itsEC;
    SpinGroup                    itsSpinGroup;     //the imposed spin subgroup (V1.37)
    qchem::Ortho                 itsBasisOrtho;    //S-orthogonalisation mode for the generalised eigenproblem
    double                       itsBasisOrthoTol; //near-null eigen/singular-value cutoff (Eigen/SVD; 0 = keep all)
    ReservoirPartition           itsPartition;     //how the EC pools its electrons (V1.11 increment 2)
    SCFAccelerator*              itsAccelerator;   // NON-template manager (RealComplexPlan §6)
    EnergyLevels                 itsELevels;
    std::map<Spin,EnergyLevels>  itsSpin_ELevels;
    std::map<Spin,std::map<Irrep,double>> itsAufbauNe; //per-irrep electron count, keyed by irrep (recomputed each iteration)

    std::vector<iwf_child_t>                itsIWFs;   // the §4 child slot: per-block scalar
    std::map<Irrep,iwf_ref_t>               itsQNWFs;  //sort by Irrep for easy lookup.
    std::map<Spin,std::vector<iwf_ref_t>>   itsSpinWFs; //Sort by spin.
};

using CompositeWF  = tCompositeWF<double>;
using cCompositeWF = tCompositeWF<dcmplx>;

} //namespace

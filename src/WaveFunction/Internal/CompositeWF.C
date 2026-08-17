// File: CompositeWF.H  Wave function as a list of Irrep wave functions.
module;
#include <vector>
#include <map>
#include <memory>
export module qchem.WaveFunction.Internal.CompositeWF;
export import qchem.WaveFunction.SCF;
import qchem.SCFAccelerator;
import qchem.WaveFunction.Internal.IrrepWF;
import qchem.WaveFunction.Types;
import qchem.LASolver;   // qchem::Ortho (the basis-overlap orthogonalisation: Cholesky | Eigen | SVD + tol)
import qchem.Matrix3D;   // Matrix3D -- the reciprocal point ops passed to each composite density (IBZ symmetrization)

export namespace qchem::WaveFunction
{

using SCFAccelerators::tSCFAccelerator;
using qchem::Hamiltonian::tHamiltonian;
using ChargeDensity::tDM_CD;
using ChargeDensity::tChargeDensity;

// Wave function as a list of per-irrep wave functions.  Templated on the matrix element type T
// (rX/cX); CompositeWF is the <double> alias (atoms/molecules), cCompositeWF the <dcmplx>
// (plane-wave / single-k Bloch-irrep) instantiation.
template <class T> class tCompositeWF
    : public virtual tSCFWaveFunction<T>
{
public:
    typedef typename tWaveFunction<T>::iqns_t iqns_t;
    using tWaveFunction<T>::GetChargeDensity;   // keep the no-arg (whole-density) overload visible past GetChargeDensity(Spin)

    //! \a basisOrtho selects how the (per-irrep) orbital-overlap S is orthogonalised for the generalised
    //! eigenproblem: \c Cholesky (default; requires S positive-definite) or \c Eigen / \c SVD with a
    //! \a basisOrthoTol cutoff that DROPS near-null eigen/singular values -- canonical orthogonalisation for a
    //! linearly-dependent basis (e.g. diffuse Gaussians on a dense lattice).  \a basisOrthoTol\f$\le0\f$ = keep all.
    tCompositeWF(const tbs_t<T>*,const ElectronConfiguration*,tSCFAccelerator<T>*,
                 qchem::Ortho basisOrtho=qchem::Auto, double basisOrthoTol=0.0);
    ~tCompositeWF();

    virtual void            DoSCFIteration  (tHamiltonian<T>&,const tChargeDensity<T>*   )      ;
    //! Iteration-0 seed.  BUILDS and hands over the first real density (V1.25: owning return).
    virtual std::unique_ptr<tDM_CD<T>> Init(tHamiltonian<T>&,const tChargeDensity<T>*,
                                            OccupationPolicy<T>&, double mergeTol);
    virtual bool            BuildFockAndComputeSteps(tHamiltonian<T>&,const tChargeDensity<T>*);
    virtual void            MoveOrbitals    (OccupationPolicy<T>&, double t, bool commit, double mergeTol,
                                             bool holdBlock=false);
    virtual const Orbitals* GetOrbitals     (const Irrep&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;}
    virtual void            FillOrbitals    (OccupationPolicy<T>&, double mergeTol, bool holdBlock=false);
    // (SetMOM/SetSmearing/GetEntropyTerm/AdoptMOMReference/ReleaseMOMReference are GONE -- the
    //  SCFIterator's OccupationPolicy slot owns that configuration and state, V1.11 inc 3.)
    virtual iqns_t          GetQNs          () const;
    virtual void            EmitBasisUsage  () const;

    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity(Spin) const;   //!< BUILDS it (V1.25)
    virtual EnergyLevels    GetEnergyLevels (Spin) const;


protected:
    void MakeIrrepWFs(Spin);

private:
    typedef tIrrepWF<T> iwf_t;
    typedef std::unique_ptr<iwf_t> uiwf_t;
    //! RANKED integer fill of one reservoir (the molecular cross-irrep aufbau, one spin channel): pick which
    //! orbitals across the reservoir's blocks are occupied, then fill each block with its resulting count.
    void FillReservoirRanked    (OccupationPolicy<T>&, const std::vector<iwf_t*>&, double mergeTol, bool useMOM);
    //! Smeared fill of one reservoir: solve ONE μ over the reservoir's blocks -- across the Bloch mesh when
    //! the partition spans spatial (doc/GPWPlan1.md item 3), across SPIN when it spans spin (free moment).
    void FillReservoirAtSharedMu(OccupationPolicy<T>&, const std::vector<iwf_t*>&, double mergeTol);

    const tbs_t<T>*              itsBS;
    const ElectronConfiguration* itsEC;
    qchem::Ortho                 itsBasisOrtho;    //S-orthogonalisation mode for the generalised eigenproblem
    double                       itsBasisOrthoTol; //near-null eigen/singular-value cutoff (Eigen/SVD; 0 = keep all)
    ReservoirPartition           itsPartition;     //how the EC pools its electrons (V1.11 increment 2)
    tSCFAccelerator<T>*          itsAccelerator;
    EnergyLevels                 itsELevels;
    std::map<Spin,EnergyLevels>  itsSpin_ELevels;
    std::map<Spin,std::map<Irrep,double>> itsAufbauNe; //per-irrep electron count, keyed by irrep (recomputed each iteration)

    std::vector<uiwf_t>                itsIWFs;
    std::map<Irrep,iwf_t*>             itsQNWFs; //sort by Irrep for easy lookup.
    std::map<Spin,std::vector<iwf_t*>> itsSpinWFs; //Sort by spin.
};

using CompositeWF  = tCompositeWF<double>;
using cCompositeWF = tCompositeWF<dcmplx>;

} //namespace

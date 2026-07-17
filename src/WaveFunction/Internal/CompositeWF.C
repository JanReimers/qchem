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
                 qchem::Ortho basisOrtho=qchem::Cholesky, double basisOrthoTol=0.0);
    ~tCompositeWF();

    virtual void            DoSCFIteration  (tHamiltonian<T>&,const tChargeDensity<T>*   )      ;
    virtual tDM_CD<T>*      Init            (tHamiltonian<T>&,const tChargeDensity<T>*, double mergeTol);
    virtual bool            BuildFockAndComputeSteps(tHamiltonian<T>&,const tChargeDensity<T>*);
    virtual void            MoveOrbitals    (double t, bool commit, double mergeTol);
    virtual const Orbitals* GetOrbitals     (const Irrep&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;}
    virtual void            FillOrbitals    (double mergeTol);
    virtual void            SetMOM          (bool useMOM, int startIter);
    virtual iqns_t          GetQNs          () const;

    virtual tDM_CD<T>*      GetChargeDensity(Spin) const;
    virtual EnergyLevels    GetEnergyLevels (Spin) const;


protected:
    void MakeIrrepWFs(Spin);

private:
    typedef tIrrepWF<T> iwf_t;
    typedef std::unique_ptr<iwf_t> uiwf_t;
    void FillOrbitalsAufbau(double mergeTol); //fill globally-lowest orbitals across all irreps

    const tbs_t<T>*              itsBS;
    const ElectronConfiguration* itsEC;
    qchem::Ortho                 itsBasisOrtho;    //S-orthogonalisation mode for the generalised eigenproblem
    double                       itsBasisOrthoTol; //near-null eigen/singular-value cutoff (Eigen/SVD; 0 = keep all)
    bool                         itsAufbau;   //molecular aufbau across irreps (vs fixed per-irrep EC)
    bool                         itsUseMOM=false;    //Maximum Overlap Method for this run (from SCFParams::UseMOM)
    bool                         itsMOMActive=false; //cross-irrep MOM armed (parked molecular path; set after 1st fill)
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

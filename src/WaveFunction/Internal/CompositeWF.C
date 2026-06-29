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

    tCompositeWF(const tbs_t<T>*,const ElectronConfiguration*,tSCFAccelerator<T>*);
    ~tCompositeWF();

    virtual void            DoSCFIteration  (tHamiltonian<T>&,const tChargeDensity<T>*   )      ;
    virtual bool            BuildFockAndComputeSteps(tHamiltonian<T>&,const tChargeDensity<T>*);
    virtual void            MoveOrbitals    (double t, bool commit, double mergeTol);
    virtual const Orbitals* GetOrbitals     (const Irrep&) const;
    virtual       Orbitals* GetOrbitals     (const Irrep&)      ;
    virtual EnergyLevels    GetEnergyLevels () const {return itsELevels;}
    virtual void            FillOrbitals    (double mergeTol);
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
    bool                         itsAufbau;   //molecular aufbau across irreps (vs fixed per-irrep EC)
    bool                         itsMOMActive=false; //once the accelerator engages, pick occupation by MOM (overlap) not eigenvalue
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

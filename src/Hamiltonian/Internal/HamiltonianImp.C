// File: HamiltonianImp.C  General matrix implementation of a Hamiltonian operator.
module;
#include <vector>
#include <memory>
export module qchem.Hamiltonian.Internal.Hamiltonian;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Structure;

export namespace qchem::Hamiltonian
{

// Templated on T (rX/cX convention): the bare HamiltonianImp is the <double> alias the existing concrete
// Hamiltonians derive from; cHamiltonianImp is the dcmplx (plane-wave) instantiation.
template <class T> class tHamiltonianImp
    : public virtual tHamiltonian<T>
{
public:
    tHamiltonianImp();
    virtual void Add(   tStatic_HT<T>* );
    virtual void Add(  tDynamic_HT<T>*);
    virtual void Add(tDynamic_HF_HT<T>*);

    using tHamiltonian<T>::GetMatrix;   // keep the base 3-arg (null-basis) convenience overload visible
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>*,const Spin& S,const tChargeDensity<T>*,const tbs_t<T>* wholeBasis);
    virtual EnergyBreakdown GetTotalEnergy  (const tDM_CD<T>* ) const;
    virtual bool            IsPolarized() const {return itsIsPolarized;}
    virtual bool            IsRelativistic() const {return itsIsRelativistic;}
    virtual std::ostream&   Write(std::ostream&) const;

protected:
    typedef std::shared_ptr<const Structure> st_t;
    void InsertStandardTerms(const st_t & st);   // molecular standard terms (double only)
    typedef std::vector<std::unique_ptr<   tStatic_HT<T>>> shtv_t;
    typedef std::vector<std::unique_ptr<  tDynamic_HT<T>>> dhtv_t;
    typedef std::vector<std::unique_ptr<tDynamic_HF_HT<T>>> ehtv_t;

    shtv_t itsSHTs;
    dhtv_t itsDHTs;
    ehtv_t itsEHTs;   // whole-system HF (exact 4-index) terms

    bool   itsIsPolarized;
    bool   itsIsRelativistic;
};

using HamiltonianImp  = tHamiltonianImp<double>;
using cHamiltonianImp = tHamiltonianImp<dcmplx>;

} //namespace
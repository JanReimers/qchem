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
    //! A density MATRIX (not a matrix-free fit) is required to seed the SCF iff this Hamiltonian either holds
    //! an exact-exchange HF term (K needs D) or is relativistic (the non-rel LDA-sibling SAD bootstrap does
    //! not apply -- see SCFIterator).  Derived from the term lists, so no concrete Hamiltonian need declare it.
    virtual bool            RequiresDensityMatrix() const {return !itsHF_HTs.empty() || itsIsRelativistic;}
    virtual std::ostream&   Write(std::ostream&) const;

protected:
    typedef std::shared_ptr<const Structure> st_t;
    void InsertStandardTerms(const st_t & st);   // molecular standard terms (double only)
    typedef std::vector<std::unique_ptr<   tStatic_HT<T>>> shtv_t;
    typedef std::vector<std::unique_ptr<  tDynamic_HT<T>>> dhtv_t;
    typedef std::vector<std::unique_ptr<tDynamic_HF_HT<T>>> hf_htv_t;

    shtv_t   itsSHTs;
    dhtv_t   itsDHTs;
    hf_htv_t itsHF_HTs;   // HF capable terms require a widened interface for efficient J/K table handling.

    bool   itsIsPolarized;
    bool   itsIsRelativistic;
};

using HamiltonianImp  = tHamiltonianImp<double>;
using cHamiltonianImp = tHamiltonianImp<dcmplx>;

} //namespace
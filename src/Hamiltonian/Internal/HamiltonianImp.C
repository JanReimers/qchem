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

// Templated on T (rX/cX convention): the bare rHamiltonianImp is the <double> alias the existing concrete
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
    typedef std::vector<std::unique_ptr<   tStatic_HT<T>>> shtv_t;
    typedef std::vector<std::unique_ptr<  tDynamic_HT<T>>> dhtv_t;
    typedef std::vector<std::unique_ptr<tDynamic_HF_HT<T>>> hf_htv_t;

    shtv_t   itsSHTs;
    dhtv_t   itsDHTs;
    hf_htv_t itsHF_HTs;   // HF capable terms require a widened interface for efficient J/K table handling.

    bool   itsIsPolarized;
    bool   itsIsRelativistic;
};

//! \brief The MOLECULAR (real) Hamiltonian implementation: \c tHamiltonianImp<double> plus the standard
//! molecular term set.  A CLASS, not an alias, so that \c InsertStandardTerms exists ONLY on the lineage
//! that has one (R2.8).
//!
//! It used to be a member of the T-generic base, which forced a \c dcmplx specialization whose entire body
//! was \c assert(false) -- the plane-wave Hamiltonian assembles its terms explicitly and has no "standard"
//! set (no bare \c Ven: the nuclei are pseudised; no molecular \c Kinetic instantiation).  Declaring it here
//! makes that a COMPILE error rather than a runtime assert, and removes an LSP-violating stub from the base.
class rHamiltonianImp : public tHamiltonianImp<double>
{
protected:
    //! Kinetic + ion-ion + bare nuclear attraction -- the terms every molecular/atomic Hamiltonian starts with.
    void InsertStandardTerms(const st_t& st);
};

using cHamiltonianImp = tHamiltonianImp<dcmplx>;

} //namespace
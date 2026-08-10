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
    //! CONJUNCTIVE over the terms -- see the base declaration.  One non-Coulombic term (a PP) invalidates
    //! the virial for the whole Hamiltonian, so this is AND-ed in Add() where the two flags above are OR-ed.
    virtual bool            IsVirialValid () const {return itsIsVirialValid;}
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
    bool   itsIsVirialValid=true;   //!< AND of the terms' IsVirialValid() (see Add); true for an empty H
};

// NO InsertStandardTerms.  There is no "standard" term set to insert (R2.8; user 2026-08-07: the whole
// idea was too clever by half).  Every concrete Hamiltonian now spells its own core out in three Adds,
// because the core genuinely varies along axes the one helper could not express:
//   * BARE vs PSEUDISED nuclei decides BOTH the ion charge (Z vs Zion) AND the electron-nuclear term(s) --
//     and that is NOT the double/dcmplx axis it was filed under: Ham_PP is <double> and never called the
//     helper either.
//   * the BASIS lineage then decides WHICH electron-nuclear implementation (PP_Local's mesh quadrature vs
//     Ven_PP_Short/_Long's G-space route) -- a second, independent axis.
//   * the Dirac Hamiltonians replace the kinetic term outright (DiracKinetic + RestMass) and carry no
//     ion-ion at all.
// Seven of eleven Hamiltonians shared the helper; for the other four it was noise, and its <dcmplx>
// specialization was a bare assert(false).  Three explicit Adds read better than a name that hides them.
using rHamiltonianImp = tHamiltonianImp<double>;
using cHamiltonianImp = tHamiltonianImp<dcmplx>;

} //namespace
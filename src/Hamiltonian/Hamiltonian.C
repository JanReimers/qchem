// File: Hamiltonian.C  Interface a Hamiltonianian operator.
export module qchem.Hamiltonian;
export import qchem.ChargeDensity;
import qchem.Streamable;
export import qchem.Energy;
export import qchem.Hamiltonian.Types;


export namespace qchem::Hamiltonian
{

using ChargeDensity::tStatic_CC;
using ChargeDensity::tDynamic_CC;
using ChargeDensity::tChargeDensity;
using ChargeDensity::rChargeDensity;
using ChargeDensity::cChargeDensity;
using ChargeDensity::tDM_CD;
using ChargeDensity::rDM_CD;
using ChargeDensity::cDM_CD;
using ChargeDensity::DM_CD;

//
//  Abstract base for any HamiltonianTerm (HT) terms in the Hamiltonian.
//  We have two distinct types of HT:
//  Static_HT - Does not depend on the Charge Dnesity (CD), and therefore does change during iterations
//  Dynamic_HT - Vee, Vac which depend on the CD, and change with each iteration.
//
//  Templated on the matrix element type T (double for atoms/molecules; dcmplx for the plane-wave
//  lattice lineage).  hmat_t<double> IS rsmat_t and tobs_t<double> IS obs_t, so the <double> aliases
//  below leave all existing real code unchanged.
//
template <class T> class tStatic_HT
    : public virtual Streamable
    , public virtual tStatic_CC<T>
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&) const=0;
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    virtual bool             IsPolarized   () const {return false;}
    virtual bool             IsRelativistic() const {return false;}
    virtual bool             RequiresDensityMatrix() const {return false;}
};

template <class T> class tDynamic_HT
    : public virtual Streamable
    , public virtual tDynamic_CC<T>
{
public:
    // Fock build consumes the DFT (matrix-free) density face; energy needs the density matrix.
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*) const=0;
    //! Cross-irrep-aware Fock build.  Default: ignore the context and defer to the one-irrep form above --
    //! so every term is a valid consumer for free, and only a term that genuinely exploits the global view
    //! (the Coulomb/exchange term, banking ERI4 bra-ket symmetry -- doc/ERI4Rework.md \S5.4) overrides it.
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s,const tChargeDensity<T>* cd,
                                       const tHamiltonianContext<T>&) const
    { return GetMatrix(bs,s,cd); }
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    virtual bool             IsPolarized   () const {return false;}
    virtual bool             IsRelativistic() const {return false;}
    virtual bool             RequiresDensityMatrix() const {return false;}
};

template <class T> class tHamiltonian
    : public virtual Streamable
{
public:
    virtual void            Add             ( tStatic_HT<T>*)=0;
    virtual void            Add             (tDynamic_HT<T>*)=0;
    //! Assemble the Fock/Hamiltonian for one irrep \a bs, given the cross-irrep \a ctx (threaded to the
    //! dynamic terms).  This is the primary form the SCF (CompositeWF/IrrepWF) drives.
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>*,const tHamiltonianContext<T>&)=0;
    //! Convenience for callers with no cross-irrep view (e.g. stand-alone tests): builds with an empty
    //! context, so every dynamic term takes its default (context-ignoring) path.
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>* bs,const Spin& s,const tChargeDensity<T>* cd)
    { return GetMatrix(bs,s,cd,tHamiltonianContext<T>{}); }
    virtual EnergyBreakdown GetTotalEnergy  (  const tDM_CD<T>*    ) const=0;
    virtual bool            IsPolarized   () const=0;
    virtual bool            IsRelativistic() const=0;
    //! DFT/KS: the Fock is a functional of rho(r) alone -> false (can be seeded from a numeric ScalarFunction).
    //! HF/DHF override -> true: they need the density MATRIX D for exact exchange K, so the SCFIterator must
    //! bootstrap them (route rho through a DFT sibling to manufacture a D0).  See project_numericcd_refactor.
    virtual bool            RequiresDensityMatrix() const {return false;}
};

// r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t); bare names transitional (= r*), rename pinned.
using rStatic_HT   = tStatic_HT<double>;   using cStatic_HT   = tStatic_HT<dcmplx>;
using rDynamic_HT  = tDynamic_HT<double>;  using cDynamic_HT  = tDynamic_HT<dcmplx>;
using rHamiltonian = tHamiltonian<double>; using cHamiltonian = tHamiltonian<dcmplx>;
using Static_HT    = rStatic_HT;
using Dynamic_HT   = rDynamic_HT;
using Hamiltonian  = rHamiltonian;

} //namespace


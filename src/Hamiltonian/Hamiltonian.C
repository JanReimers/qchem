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
};

template <class T> class tDynamic_HT
    : public virtual Streamable
    , public virtual tDynamic_CC<T>
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>*,const Spin&,const tDM_CD<T>*) const=0;
    virtual void             GetEnergy(EnergyBreakdown&,  const tDM_CD<T>*) const=0;
    virtual bool             IsPolarized   () const {return false;}
    virtual bool             IsRelativistic() const {return false;}
};

template <class T> class tHamiltonian
    : public virtual Streamable
{
public:
    virtual void            Add             ( tStatic_HT<T>*)=0;
    virtual void            Add             (tDynamic_HT<T>*)=0;
    virtual hmat_t<T>       GetMatrix(const tobs_t<T>*,const Spin&,const tDM_CD<T>*)=0;
    virtual EnergyBreakdown GetTotalEnergy  (  const tDM_CD<T>*    ) const=0;
    virtual bool            IsPolarized   () const=0;
    virtual bool            IsRelativistic() const=0;
};

// r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t); bare names transitional (= r*), rename pinned.
using rStatic_HT   = tStatic_HT<double>;   using cStatic_HT   = tStatic_HT<dcmplx>;
using rDynamic_HT  = tDynamic_HT<double>;  using cDynamic_HT  = tDynamic_HT<dcmplx>;
using rHamiltonian = tHamiltonian<double>; using cHamiltonian = tHamiltonian<dcmplx>;
using Static_HT    = rStatic_HT;
using Dynamic_HT   = rDynamic_HT;
using Hamiltonian  = rHamiltonian;

} //namespace


// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
#include <cassert>
export module qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Orbital_1E_IBS;
import qchem.BasisSet1.DB_Cache;
import qchem.BasisSet.Internal.ERI4;

export namespace BasisSet1
{

template <class T> class Orbital_HF_IBS;
//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
template <class T> class Integrals_HF : public virtual IrrepBasisSet_IDs
{
public:
    public:
    virtual ERI4  MakeDirect  (const Orbital_HF_IBS<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    virtual ERI4  MakeExchange(const Orbital_HF_IBS<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    const   ERI4& Direct(const Orbital_HF_IBS<T>& c) const
    {
        auto cache=theGlobalCache;
        assert(cache);
        IntegralsCache_Base::IBS_ID_t ab(RadialID(),AngularID());
        IntegralsCache_Base::IBS_ID_t cd(c.RadialID(),c.AngularID());
        return cache->Has(IntegralsCache_Base::I4C::Direct,ab,cd)
            ? cache->GetERI4() : cache->SetDirect(MakeDirect(c));
    }
    const   ERI4& Exchange(const Orbital_HF_IBS<T>& c) const
    {
        auto cache=theGlobalCache;
        assert(cache);
        IntegralsCache_Base::IBS_ID_t ab(RadialID(),AngularID());
        IntegralsCache_Base::IBS_ID_t cd(c.RadialID(),c.AngularID());
        return cache->Has(IntegralsCache_Base::I4C::Exchange,ab,cd)
            ? cache->GetERI4() : cache->SetExchange(MakeExchange(c));
    }


};

//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs).
//! 
template <class T> class Integrals_BS_HF
{
public:
    virtual ERI4 Direct  (const Orbital_HF_IBS<T>& a,const Orbital_HF_IBS<T>& c) const=0;
    virtual ERI4 Exchange(const Orbital_HF_IBS<T>& a,const Orbital_HF_IBS<T>& b) const=0;
};


template <class T> class Orbital_HF_IBS
    : public virtual Integrals_HF<T> //Two electron integrals used for HF
{
public:
    using Integrals_HF<T>::Direct;
    using Integrals_HF<T>::Exchange;
    virtual void AccumulateDirect  (rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
    virtual void AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
};

} //namespace

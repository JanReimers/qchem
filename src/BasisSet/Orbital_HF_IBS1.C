// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
#include <cassert>
export module qchem.Orbital_HF_IBS1;
export import qchem.Orbital_1E_IBS1;
import qchem.BasisSet.DB_Cache1;
import qchem.BasisSet.Internal.ERI4;

export template <class T> class Orbital_HF_IBS1;
//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
export template <class T> class Integrals_HF1 : public virtual IrrepBasisSet_IDs
{
public:
    public:
    virtual ERI4  MakeDirect  (const Orbital_HF_IBS1<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    virtual ERI4  MakeExchange(const Orbital_HF_IBS1<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    const   ERI4& Direct(const Orbital_HF_IBS1<T>& c) const
    {
        auto cache=theGlobalCache;
        assert(cache);
        IntegralsCache_Base::IBS_ID_t ab(RadialID(),AngularID());
        IntegralsCache_Base::IBS_ID_t cd(c.RadialID(),c.AngularID());
        return cache->Has(IntegralsCache_Base::I4C::Direct,ab,cd)
            ? cache->GetERI4() : cache->SetDirect(MakeDirect(c));
    }
    const   ERI4& Exchange(const Orbital_HF_IBS1<T>& c) const
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
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs)
//! using IBS IDs.
// export template <class T> class Integrals_BS_HF
// {
// public:
//     typedef UniqueID::IDtype IDType;   
//     virtual ERI4 Direct  (IDType a,IDType c) const=0;
//     virtual ERI4 Exchange(IDType a,IDType b) const=0;
// };


export template <class T> class Orbital_HF_IBS1
    : public virtual Orbital_1E_IBS1<T>
    , public virtual Integrals_HF1<T> //Two electron integrals used for HF
{
public:
    Orbital_HF_IBS1(const Irrep_QNs::sym_t& sym) : Orbital_1E_IBS1<T>(sym) {};
    using Integrals_HF1<T>::Direct;
    using Integrals_HF1<T>::Exchange;
    virtual void AccumulateDirect  (rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS1<T>* bs_cd) const;
    virtual void AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS1<T>* bs_cd) const;

};


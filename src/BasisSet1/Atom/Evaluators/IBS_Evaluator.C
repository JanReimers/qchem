// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <vector>
#include <ranges>
#include <iosfwd>
#include <cassert>
#include <sstream>
#include "forward.H"
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.Evaluators.IBS;
export import qchem.BasisSet1.Atom.Evaluators.Internal.ExponentGrouper;
import qchem.BasisSet1.Internal.Cache4;
import qchem.BasisSet1.Atom.Evaluators.Internal.AngularIntegrals;
export import qchem.BasisSet1.Internal.ERI3;
export import qchem.Symmetry.Irrep;
export import qchem.VectorFunction;
import qchem.Symmetry.Ylm;

export using dERI3=ERI3<double>;

export template <class E> concept isGeneric_Evaluator = requires (const E& e,size_t i, size_t j, const rvec3_t& r)
{
    e.size();
    e.operator()(r);
    e.Gradient  (r);
    e.Norm     (i);
};

export template <class E> concept is1E_Evaluator = isGeneric_Evaluator<E> && requires  (E e,size_t i, size_t j, const rvec3_t& r)
{
    e.Norm     (i);
    e.Overlap(i,j); //Should all be inline.
    e.Grad2  (i,j);
    e.Inv_r1 (i,j);
    e.Inv_r2 (i,j);
};

export template <class E> concept isFit_Evaluator = isGeneric_Evaluator<E> && requires  (E e,size_t i, size_t j, size_t ic)
{
    e.Norm     (i);
    e.Overlap(i,j); //Should all be inline.
    e.Repulsion(i,j);
    e.Charge   (i);
};
// Support 3C Overlap and Repulsion
export template <class E> concept isDFT_Evaluator = requires (E e,size_t i, size_t j, size_t ic)
{
    e.Overlap  (i,j,e,ic); 
    e.Repulsion(i,j,e,ic);
};
// Support 4C Hartree-Fock Direct/Exchange integrals.
export template <class E> concept isHF_Evaluator = isGeneric_Evaluator<E> && requires (E a,size_t l,const Cacheable* c ,AngularIntegrals::rvec11_t Ak)
{
    a.maxSpan();
    // a.size();
    a.Getl();
    a.RadialType(); 
    a.indices();
    a.MakeCache4();
    a.direct(c,l,l,Ak);
    a.exchange(c,l,l,Ak);
};

// Support cross Kinetic
export template <class E> concept isRKBL_Evaluator = is1E_Evaluator<E> && requires  (E e,size_t i, size_t j)
{
    e.Grad2    (i,j,e);
    e.Inv_r2   (i,j,e);
};


export template <class E> concept isFull_NR_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E> && isDFT_Evaluator<E>;
export template <class E> concept isHF_NR_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E>;
export template <class E> concept isFull_HF_Evaluator = isGeneric_Evaluator<E> && is1E_Evaluator<E> && isDFT_Evaluator<E> && isHF_NR_Evaluator<E>;


export class IBS_Evaluator 
    : public virtual Cache4_Client
    , public VectorFunction<double>
{
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
public:
    using is_t=std::vector<int>;
    

    IBS_Evaluator(int _l, const is_t& _mls) : l(_l), mls(_mls),ns(0),grouper(0) {};
    IBS_Evaluator(const Irrep_QNs::sym_t& ylm);
    virtual ~IBS_Evaluator() {};

    virtual void          Register(Grouper*)=0; //Set up unique spline or exponent indexes.
    virtual size_t        size    (             ) const {return ns.size();}
    virtual int           Getl    (             ) const {return l;}
    virtual size_t        es_index(size_t i     ) const {return es_indices[i];}
    virtual const is_t&   Getmls  (             ) const {return mls;}

    iota_view             indices (             ) const {return iota_view(size_t(0),size());}
    iota_view             indices (size_t start ) const {return iota_view(start,size());}
    virtual std::ostream& Write   (std::ostream&) const=0;
    virtual size_t maxSpan() const {return size();}  //assume no overlap for indeces separated by > maxSpan

    virtual size_t        GetVectorSize() const {return size();}

    virtual  rvec_t Norm      () const=0;
   
    virtual std::string RadialID () const=0;
    virtual std::string AngularID() const;
    virtual std::string Name() const=0;
    static AngularIntegrals::rvec11_t Coulomb_AngularIntegrals(const IBS_Evaluator& a,const IBS_Evaluator& c)
    {
        AngularIntegrals::rvec11_t Ak(0.0);
        int la=a.Getl(),lc=c.Getl();
        const IBS_Evaluator::is_t& amls=a.Getmls(),cmls=c.Getmls();
        size_t nac=amls.size()*cmls.size();
        size_t g=(2*la+1)*(2*lc+1); //degenracy
        if (nac==g || nac==0)
        {
            Ak+=AngularIntegrals::Coulomb(la,lc);
        }
        else
        {
            // For direct integrals these actually factor.  But for exchange they do not.
            // So it may not be worth while coding the factored version here.
            for (auto ma:amls)
            for (auto mc:cmls)
                Ak+=AngularIntegrals::Coulomb(la,lc,ma,mc);
            Ak/=(double)nac;
        }
        return Ak;
    }
    static AngularIntegrals::rvec11_t ExchangeAngularIntegrals(const IBS_Evaluator& a,const IBS_Evaluator& b)
    {
        AngularIntegrals::rvec11_t Ak(0.0);
        int la=a.Getl(),lb=b.Getl();
        const IBS_Evaluator::is_t& amls=a.Getmls(),bmls=b.Getmls();
        size_t nab=amls.size()*bmls.size();
        size_t g=(2*la+1)*(2*lb+1); //degenracy
        if (nab==0 || nab==g)
        {
            Ak+=AngularIntegrals::Exchange(la,lb);
        }
        else
        {
            for (auto ma:amls)
            for (auto mb:bmls)
                Ak+=AngularIntegrals::Exchange(la,lb,ma,mb);
            Ak/=nab;
        }
        return Ak;
    }    

protected:
    friend class Cache4Tests;

    int  l;
private:
    // used by Coulomb_AngularIntegrals and ExchangeAngularIntegrals
    is_t   mls;
protected:
    rvec_t ns;
    const ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index

};



IBS_Evaluator::IBS_Evaluator(const Irrep_QNs::sym_t& y) :  l(0), mls({}),ns(0),grouper(0)
{
    const Yl_Sym* yl=dynamic_cast<const Yl_Sym*>(y.get());
    assert(yl);
    l=yl->GetL();
    const Ylm_Sym* ylm=dynamic_cast<const Ylm_Sym*>(y.get());
    if (ylm)
    {
        mls=ylm->Getmls();
    }
}

std::string IBS_Evaluator::AngularID() const
{
     std::ostringstream os;
     os << l << " {";
     for (auto ml:mls) os << ml << " ";
     os << "}";
     return os.str();
}
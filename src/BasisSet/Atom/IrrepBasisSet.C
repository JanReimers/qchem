// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include <cassert>
#include <iostream>
#include "forward.H"

export module qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Atom.Evaluators;
import qchem.BasisSet.Internal.Cache4;
import qchem.Symmetry.Factory;
import qchem.Constants;
import qchem.Blaze;

export namespace qchem::BasisSet
{
namespace Atom
{

    using namespace Evaluators;
    using EvaluatorsBase = Evaluators::Evaluator;
//
//  Common IrrepBasisSet functionality for atom basis sets.  All the work is done by the evaluator
//
template <isOpr_Evaluator E> class IrrepBasisSetImp
    : public virtual BasisSet::IrrepBasisSet<double>
    , public virtual IrrepBasisSet_IDs
    , public BasisSet::IrrepBasisSetImp<double> //Pulls in Symmetry support
{
public:
    IrrepBasisSetImp(const sym_t& yl) : BasisSet::IrrepBasisSetImp<double>(yl) {}

    virtual size_t GetNumFunctions() const {return Cast().size();}
    // using statements in the final class don't seem to work, so we need to function forward.
    virtual rvec_t     operator() (const rvec3_t& r) const {return Cast().operator()(r);}
    virtual rvec3vec_t Gradient   (const rvec3_t& r) const {return Cast().Gradient  (r);}
   
    virtual std::string RadialID () const {return Cast().RadialID();}
    virtual std::string AngularID() const {return Cast().AngularID();}
    virtual std::string Name     () const {return Cast().Name();}
protected:
    auto& Cast() const {return dynamic_cast<const E&>(*this);}
};

//
//  Implement these integral engines separately as they shared between NR and RKB 1E orbital IBS implementations.
//  Fit_IBS also uses Overlap.
//
template <is1E_Evaluator E> class Integrals_Overlap
: public virtual BasisSet::Integrals_Overlap<double>
{
protected:
    virtual smat_t<double> MakeOverlap() const 
    {
        auto& e=dynamic_cast<const E&>(*this);
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Overlap(i,j);

        return S;
    }
};
template <is1E_Evaluator E> class Integrals_Kinetic
: public virtual BasisSet::Integrals_Kinetic<double>
{
protected:
    // Builds the kinetic BUILDING BLOCK \f$\langle p^2\rangle=\langle -\nabla^2\rangle\f$ (NOT energy,
    // no 1/2 -- see BasisSet/Orbital_1E_IBS.C).  For atoms \f$-\nabla^2\f$ splits into the radial
    // gradient piece e.Grad2(i,j) plus the centrifugal (angular) term \f$l(l+1)\langle r^{-2}\rangle\f$.
    virtual smat_t<double> MakeKinetic() const
    {
        auto& e=dynamic_cast<const E&>(*this);
        int l=e.Getl();
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Grad2(i,j) + l*(l+1)*e.Inv_r2(i,j);

        return S;
    }
};
template <is1E_Evaluator E> class Integrals_Nuclear
: public virtual BasisSet::Integrals_Nuclear<double>
{
protected:
    virtual smat_t<double> MakeNuclear(const Structure* cl) const 
    {
        assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
        auto& e=dynamic_cast<const E&>(*this);
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= Z*e.Inv_r1(i,j);

        return S;
    }
};

// Use E prefix to avoid name clash with the interface class Fit_IBS
template <isFit_Evaluator Evaluator> class EFit_IBS
    : public virtual Fit_IBS 
    , public Integrals_Overlap<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
    using IrrepBasisSetImp<Evaluator>::Cast;
public:
    EFit_IBS(const Evaluator& e) : IrrepBasisSetImp<Evaluator>(Symmetry::YFactory()), Evaluator(e) {};

    virtual rsmat_t MakeRepulsion(                ) const 
    {
        auto& e=Cast();
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Repulsion(i,j);

        return S;
    }
    virtual  rmat_t MakeRepulsion(const FIT_CD_ABS& f) const
    {
        auto& ea=Cast();
        auto& eb=dynamic_cast<const Evaluator&>(f);
        rmat_t S(ea.size(),eb.size());
        for (auto i:ea.indices())
            for (auto j:eb.indices())
                S(i,j)= ea.Repulsion(i,j,eb);

        return S;
    }
    virtual  rvec_t MakeCharge   (                ) const 
    {
        auto& e=Cast();
        rvec_t c(e.size());
        for (auto i:e.indices())
            c[i]=e.Charge(i);
        return c;
    }
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << "Atom fit IBS " << Name() << " " << GetSymmetry();
    }

};

//
//  1E NR orbital for atoms.  Use mixins to get the integral evaluations.
//
template <is1E_Evaluator E> class Orbital_1E_IBS
    : public virtual BasisSet::Orbital_1E_IBS<double> //This part has the symmetry.
    , public Integrals_Overlap<E>
    , public Integrals_Kinetic<E>
    , public Integrals_Nuclear<E>
{
public:
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << Name() << " " << GetSymmetry();
    }
};


template <isDFT_Evaluator E> class Orbital_DFT_IBS
    : public virtual BasisSet::Orbital_DFT_IBS<double>
{
protected:
    virtual ERI3<double> MakeOverlap3C  (const FIT_SF_ABS& _c) const
    {
        auto& ab=dynamic_cast<const E&>(*this);
        auto& c =dynamic_cast<const E&>(_c);
        ERI3<double> S3;
        for (auto ic:c.indices())
        {
            rsmat_t S(ab.size());
            for (auto i:ab.indices())
                for (auto j:ab.indices(i))
                    S(i,j)=ab.Overlap(i,j,c,ic);
            
            S3.push_back(S);
        }
        return S3;

    }
    virtual ERI3<double> MakeRepulsion3C(const FIT_CD_ABS& _c) const
    {
        auto& ab=dynamic_cast<const E&>(*this);
        auto& c =dynamic_cast<const E&>(_c);
        ERI3<double> S3;
        for (auto ic:c.indices()) 
        {
            rsmat_t S(ab.size());
            for (auto i:ab.indices())
                for (auto j:ab.indices(i))
                    S(i,j)=ab.Repulsion(i,j,c,ic);  
            
            S3.push_back(S);
        }
        return S3;

    }
};




template <isHF_Evaluator E> class Orbital_HF_IBS
    : public virtual BasisSet::Orbital_HF_IBS<double> 

{
protected:
    virtual ERI4 MakeDirect  (const BasisSet::Orbital_HF_IBS<double>& _c) const;
    virtual ERI4 MakeExchange(const BasisSet::Orbital_HF_IBS<double>& _c) const;

};

template <is1E_Evaluator E> class Orbital_RKBL_IBS
    : public virtual BasisSet::Orbital_RKBL_IBS<double> 
    , public Integrals_Overlap<E>
    , public Integrals_Kinetic<E> //Used internally only.
    , public Integrals_Nuclear<E>
{
public:
    // RELATIVISTIC kinetic (RKB L/S cross term).  Integrals_Kinetic<E>::MakeKinetic() is the
    // \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block (no 1/2); here it is scaled by 1/(2 c_light).
    // WARNING: this 1/(2 c_light) and the c_light in Imp/Orbital_DHF_IBS.C::MakeKinetic() multiply out
    // to ~1/2 -- the c's cancel.  Whether each factor sits in the *right* place is UNVERIFIED (see the
    // \warning in BasisSet/Internal/Orbital_DHF_IBS.C); kept as-is, guarded by A_DHF.
    virtual rmat_t  MakeKinetic(const Orbital_RKBS_IBS<double>& rkbs) const
    {
        auto& ea=dynamic_cast<const E&>(*this);
        auto& eb=dynamic_cast<const E&>(rkbs); //The RKBS evaluator should be derived from the RKBL evaluator.
        assert(ea.Getl()==eb.Getl());
        assert(ea.size()==eb.size());
        return Integrals_Kinetic<E>::MakeKinetic()/(2*c_light);
    }
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << Name() << " " << GetSymmetry();
    }
};

template <is1E_Evaluator E> class Orbital_RKBS_IBS
    : public virtual BasisSet::Orbital_RKBS_IBS<double> 
    , public Integrals_Kinetic<E>
    , public Integrals_Nuclear<E> //RKBS Evaluator overrides Inv_r1 definition
{
public:
    // RKB small-component OVERLAP.  In RKB the small-component metric \f$\langle\chi^S|\chi^S\rangle\f$
    // is proportional to \f$\langle p^2\rangle\f$, so this reuses the kinetic \f$\langle p^2\rangle\f$
    // block (Integrals_Kinetic::MakeKinetic, no 1/2) as the overlap.  Any \f$1/(4c^2)\f$ normally
    // expected here is absent -- presumed to cancel against RestMass/scaling elsewhere; UNVERIFIED
    // (see the \warning in BasisSet/Internal/Orbital_DHF_IBS.C).
    virtual rsmat_t MakeOverlap() const
    {
        return Integrals_Kinetic<E>::MakeKinetic();
    }
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << Name() << " " << GetSymmetry();
    }

};



template <isHF_Evaluator E> ERI4 Orbital_HF_IBS<E>::MakeDirect(const BasisSet::Orbital_HF_IBS<double>& _c) const 
{
    auto& a=dynamic_cast<const E&>(*this);
    auto& c=dynamic_cast<const E&>(_c);
    int la=a.Getl(), lc=c.Getl();
    assert(a.RadialType()==c.RadialType());
    size_t spanab=a.maxSpan(),spancd=c.maxSpan();
    size_t Na=a.size(), Nc=c.size();
    rvec11_t Akac=a.DirectAk(c);
    const Cache4* Rk_cache=BasisSet::theCache<double>().GetCache4(a.RadialType());
    ERI4 J(Na,Nc);

    for (size_t ia:a.indices())
    {
        Rk_cache->loop_1(a.es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c.indices())
        {
            Rk_cache->loop_2(c.es_index(ic));
            for (size_t ib:a.indices())
            {
                if (ib<ia) continue; 
                if (ib>ia+spanab) continue;
                rsmat_t& Jab=J(ia,ib);
                Rk_cache->loop_3(a.es_index(ib));
                for (size_t id:c.indices())
                {
                    if (id<ic) continue;
                    if (id>ic+spancd) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        std::cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        std::cout << Jab(ic,id) << std::endl;    
                        assert(false);
                    }
                    double norm=a.Norm(ia)*a.Norm(ib)*c.Norm(ic)*c.Norm(id);
                    const Cacheable4* rk=Rk_cache->loop_4(c.es_index(id));
                    double AkacRkac=E::direct(rk,la,lc,Akac);  //static call.
                    Jab(ic,id)=AkacRkac*norm;
                }
            }
        }
    }
    return J;
}

template <isHF_Evaluator E> ERI4 Orbital_HF_IBS<E>::MakeExchange(const BasisSet::Orbital_HF_IBS<double>& _c) const 
{
    auto& a=dynamic_cast<const E&>(*this);
    auto& c=dynamic_cast<const E&>(_c);
    assert(a.RadialType()==c.RadialType());
    size_t spanab=a.maxSpan(),spancd=c.maxSpan();
    size_t Na=a.size(), Nc=c.size();
    int la=a.Getl(), lc=c.Getl();
    rvec11_t Akac=a.ExchangeAk(c);
    const Cache4* Rk_cache=BasisSet::theCache<double>().GetCache4(a.RadialType());

    ERI4 K(Na,Nc);
    for (size_t ia:a.indices())
    {
        Rk_cache->loop_1(a.es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c.indices())
        {
            if (ia>ic+spanab || ic>ia+spancd) continue;
            for (size_t ib:a.indices(ia))
            {
                rsmat_t& Kab=K(ia,ib);
                Rk_cache->loop_2(a.es_index(ib));
                Rk_cache->loop_3(c.es_index(ic));
                for (size_t id:c.indices())
                {
                    if (id>ib+spancd || ib>id+spanab) continue;
                    double norm=a.Norm(ia)*a.Norm(ib)*c.Norm(ic)*c.Norm(id);
                    const Cacheable4* rk=Rk_cache->loop_4(c.es_index(id));
                    double AkacRKac=E::exchange(rk,la,lc,Akac); //static call.
                    if (ic==id)
                        Kab(ic,id)=AkacRKac*norm; 
                    else if (id<ic)
                        Kab(id,ic)+=0.5*AkacRKac*norm; 
                    else
                        Kab(ic,id)+=0.5*AkacRKac*norm; 

                }
            }
        }
    }

    return K;

}


}} //namespaces
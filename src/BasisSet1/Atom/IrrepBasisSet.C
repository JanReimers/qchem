// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include <cassert>
#include "forward.H"

export module qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Internal.IrrepBasisSetImp;
import qchem.BasisSet1.Internal.Orbital_DHF_IBS;
import qchem.BasisSet1.IrrepBasisSet;
import qchem.BasisSet1.Orbital_1E_IBS;
import qchem.BasisSet1.Orbital_DFT_IBS;
import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet1.Atom.Evaluators.BS;

export namespace BasisSet1
{
namespace Atom
{
//
//  Common IrrepBasisSet functionality for atom basis sets.  All the work is done by the evaluator
//
template <isEvaluator E> class IrrepBasisSetImp
    : public virtual BasisSet1::IrrepBasisSet<double>
    , public virtual IrrepBasisSet_IDs
    , public BasisSet1::IrrepBasisSetImp<double> //Pulls in Symmetry support
{
public:
    IrrepBasisSetImp(const Irrep_QNs::sym_t& yl) : BasisSet1::IrrepBasisSetImp<double>(yl) {}
    virtual rvec_t     operator() (const rvec3_t& r) const {return Cast().operator()(r);}
    virtual rvec3vec_t Gradient   (const rvec3_t& r) const {return Cast().Gradient(r);}
   
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
// 

template <isEvaluator E> class Integrals_Overlap
: public virtual BasisSet1::Integrals_Overlap<double>
{
protected:
    virtual smat_t<double> MakeOverlap() const 
    {
        auto& e=dynamic_cast<const E&>(*this);
        size_t N=e.size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= e.Overlap(i,j);

        return S;
    }
};
template <isEvaluator E> class Integrals_Kinetic
: public virtual BasisSet1::Integrals_Kinetic<double>
{
protected:
    virtual smat_t<double> MakeKinetic() const 
    {
        auto& e=dynamic_cast<const E&>(*this);
        size_t N=e.size();
        int l=e.Getl();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= e.Grad2(i,j) + l*(l+1)*e.Inv_r2(i,j);

        return S;
    }
};
template <isEvaluator E> class Integrals_Nuclear
: public virtual BasisSet1::Integrals_Nuclear<double>
{
protected:
    virtual smat_t<double> MakeNuclear(const Cluster* cl) const 
    {
        assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
        auto& e=dynamic_cast<const E&>(*this);
        size_t N=e.size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= Z*e.Inv_r1(i,j);

        return S;
    }
};

//
//  1E orbital for atoms.  Use mixins to get the integral evaluations.
//
template <isEvaluator E> class Orbital_1E_IBS
    : public virtual BasisSet1::Orbital_1E_IBS<double> //This part has the symmetry.
    , public Integrals_Overlap<E>
    , public Integrals_Kinetic<E>
    , public Integrals_Nuclear<E>
{
public:
    virtual std::ostream&  Write(std::ostream& os) const
    {
        os << "Orbital IBS " << Name() << " ";
        os << "Symmetry=" << GetSymmetry() << " ";
        auto& e=dynamic_cast<const E&>(*this);
        e.Write(os);
        return os;
    }
};


template <isEvaluator E> class Orbital_DFT_IBS
    : public virtual BasisSet1::Orbital_DFT_IBS<double>
{
protected:
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& _c) const
    {
        auto& ab=dynamic_cast<const E&>(*this);
        auto& c =dynamic_cast<const E&>(_c);
        ERI3<double> S3;
        size_t N=ab.size();
        for (size_t ic=0;ic<c.size();ic++) 
        {
            rsmat_t S(N);
            for (auto i:iv_t(0,N))
                for (auto j:iv_t(i,N))
                    S(i,j)=ab.Overlap(i,j,c,ic);  
            
            S3.push_back(S);
        }
        return S3;

    }
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& _c) const
    {
        auto& ab=dynamic_cast<const E&>(*this);
        auto& c =dynamic_cast<const E&>(_c);
        ERI3<double> S3;
        size_t N=ab.size();
        for (size_t ic=0;ic<c.size();ic++) 
        {
            rsmat_t S(N);
            for (auto i:iv_t(0,N))
                for (auto j:iv_t(i,N))
                    S(i,j)=ab.Repulsion(i,j,c,ic);  
            
            S3.push_back(S);
        }
        return S3;

    }
};



template <isEvaluator E> class Orbital_HF_IBS
    : public virtual BasisSet1::Orbital_HF_IBS<double> 

{
protected:
    Orbital_HF_IBS(BS_Evaluator* bse)  : itsEvaluator(bse) {assert(itsEvaluator);} 

    virtual ERI4 MakeDirect  (const BasisSet1::Orbital_HF_IBS<double>& _c) const 
    {
        auto& a=dynamic_cast<const E&>(*this);
        auto& c=dynamic_cast<const E&>(_c);
        return itsEvaluator->Direct(&a,&c);
    }
    virtual ERI4 MakeExchange(const BasisSet1::Orbital_HF_IBS<double>& _c) const 
    {
        auto& a=dynamic_cast<const E&>(*this);
        auto& c=dynamic_cast<const E&>(_c);
        return itsEvaluator->Exchange(&a,&c);
    }
private: 
    BS_Evaluator* itsEvaluator;
};

template <isEvaluator E> class Orbital_RKBL_IBS
    : public virtual BasisSet1::Orbital_RKBL_IBS<double> 
    , public Integrals_Overlap<E>
    , public Integrals_Nuclear<E>
{
public:
    virtual rmat_t  MakeKinetic(const Orbital_RKBS_IBS<double>& rkbs) const
    {
        auto& ea=dynamic_cast<const E&>(*this);
        auto& eb=dynamic_cast<const E&>(rkbs);
        assert(ea.Getl()==eb.Getl());
        size_t Na=ea.size(),Nb=eb.size();
        int l=ea.Getl();
        rmat_t S(Na,Nb);
        for (auto i:iv_t(0,Na))
            for (auto j:iv_t(0,Nb))
                S(i,j)= ea.Grad2(i,j,eb) + l*(l+1)*ea.Inv_r2(i,j,eb);

        return S;
    }
};

template <isEvaluator E> class Orbital_RKBS_IBS
    : public virtual BasisSet1::Orbital_RKBS_IBS<double> 
    , public Integrals_Kinetic<E>
    , public Integrals_Nuclear<E> //RKBS Evaluator overrides Inv_r1 definition
{
    virtual rsmat_t MakeOverlap() const
    {
        return Integrals_Kinetic<E>::MakeKinetic();
    }
};

}} //namespaces
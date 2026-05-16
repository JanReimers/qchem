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
import qchem.Symmetry.Yl;

export namespace BasisSet1
{
namespace Atom
{
//
//  Common IrrepBasisSet functionality for atom basis sets.  All the work is done by the evaluator
//
template <isGeneric_Evaluator E> class IrrepBasisSetImp
    : public virtual BasisSet1::IrrepBasisSet<double>
    , public virtual IrrepBasisSet_IDs
    , public BasisSet1::IrrepBasisSetImp<double> //Pulls in Symmetry support
{
public:
    IrrepBasisSetImp(const Irrep_QNs::sym_t& yl) : BasisSet1::IrrepBasisSetImp<double>(yl) {}

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
// 

template <is1E_Evaluator E> class Integrals_Overlap
: public virtual BasisSet1::Integrals_Overlap<double>
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
: public virtual BasisSet1::Integrals_Kinetic<double>
{
protected:
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
: public virtual BasisSet1::Integrals_Nuclear<double>
{
protected:
    virtual smat_t<double> MakeNuclear(const Cluster* cl) const 
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

template <isFit_Evaluator Evaluator> class Fit_IBS
    : public virtual BasisSet1::Fit_IBS 
    , public Integrals_Overlap<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
    using IrrepBasisSetImp<Evaluator>::Cast;
public:
    Fit_IBS(const Evaluator& e) : IrrepBasisSetImp<Evaluator>(Irrep_QNs::sym_t(new Yl_Sym(0))), Evaluator(e) {};

    virtual rsmat_t MakeRepulsion(                ) const 
    {
        auto& e=Cast();
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Repulsion(i,j);

        return S;
    }
    virtual  rmat_t MakeRepulsion(const BasisSet1::Fit_IBS& f) const 
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
        os << "Atom fit IBS ";
        Evaluator::Write(os);
        return os;
    }

};

//
//  1E orbital for atoms.  Use mixins to get the integral evaluations.
//
template <is1E_Evaluator E> class Orbital_1E_IBS
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


template <isDFT_Evaluator E> class Orbital_DFT_IBS
    : public virtual BasisSet1::Orbital_DFT_IBS<double>
{
protected:
    virtual ERI3<double> MakeOverlap3C  (const BasisSet1::Fit_IBS& _c) const
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
    virtual ERI3<double> MakeRepulsion3C(const BasisSet1::Fit_IBS& _c) const
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



template <isGeneric_Evaluator E> class Orbital_HF_IBS
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

template <isRKBL_Evaluator E> class Orbital_RKBL_IBS
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
        int l=ea.Getl();
        rmat_t S(ea.size(),eb.size());
        for (auto i:ea.indices())
            for (auto j:eb.indices())
                S(i,j)= ea.Grad2(i,j,eb) + l*(l+1)*ea.Inv_r2(i,j,eb);

        return S;
    }
};

template <is1E_Evaluator E> class Orbital_RKBS_IBS
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
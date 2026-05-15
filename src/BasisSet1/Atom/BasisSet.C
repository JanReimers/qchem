// File: BasisSet/Atom/BSpline/NR/BSpline_BS_Evaluator.C BSpline Basis Set for atoms.
module;
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.BasisSet;

export import qchem.BasisSet1;
export import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet1.Fit_IBS;

export import qchem.Symmetry.AtomEC;
export import qchem.Symmetry.Irrep;

import qchem.Symmetry.Yl;
import qchem.BasisSet1.Atom.Evaluators.BS;
import qchem.BasisSet1.Atom.Evaluators.IBS;
import qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Internal.BasisSetImp;
import qchem.BasisSet1.Internal.Orbital_DHF_IBS;
import qchem.BasisSet1.Internal.IrrepBasisSetImp;

export 
namespace BasisSet1 {
namespace Atom {

template <class Evaluator> class Fit_IBS
    : public virtual BasisSet1::Fit_IBS 
    , public Integrals_Overlap<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
    using IrrepBasisSetImp<Evaluator>::Cast;
public:
    Fit_IBS(const Evaluator& e) : IrrepBasisSetImp<Evaluator>(Irrep_QNs::sym_t(new Yl_Sym(0))), Evaluator(e) {};

    virtual size_t GetNumFunctions() const {return Evaluator::size();}
    virtual rsmat_t MakeRepulsion(                ) const 
    {
        auto& e=Cast();
        size_t N=e.size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= e.Repulsion(i,j);

        return S;
    }
    virtual  rmat_t MakeRepulsion(const BasisSet1::Fit_IBS& f) const 
    {
        auto& ea=Cast();
        auto& eb=dynamic_cast<const Evaluator&>(f);
        size_t Na=ea.size(),Nb=eb.size();
        rmat_t S(Na,Nb);
        for (auto i:iv_t(0,Na))
            for (auto j:iv_t(0,Nb))
                S(i,j)= ea.Repulsion(i,j,eb);

        return S;
    }
    virtual  rvec_t MakeCharge   (                ) const 
    {
        auto& e=Cast();
        size_t N=e.size();
        rvec_t c(N);
        for (auto i:iv_t(0,N))
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


template <class Evaluator> class Orbital_IBS 
    : public Orbital_1E_IBS<Evaluator>
    , public Orbital_DFT_IBS<Evaluator>
    , public Orbital_HF_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    Orbital_IBS(BS_Evaluator* bse,size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_HF_IBS<Evaluator>(bse)
    , IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,yl)
    {};
    Orbital_IBS(BS_Evaluator* bse,const rvec_t& es, const Irrep_QNs::sym_t& yl)
    : Orbital_HF_IBS<Evaluator>(bse)
    , IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(es,yl)
    {};

    virtual BasisSet1::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const 
    {
        return new Fit_IBS(Evaluator::Rescale(2.0));
    }
    virtual BasisSet1::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const
    {
        return new Fit_IBS(Evaluator::Rescale(2.0/3.0));
    }

    virtual size_t GetNumFunctions() const {return Evaluator::size();}
};

template <class Evaluator> class EOrbital_RKBL_IBS 
    : public Orbital_RKBL_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    EOrbital_RKBL_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,yl)
    {};

    virtual size_t GetNumFunctions() const {return Evaluator::size();}

    virtual std::ostream& Write(std::ostream& os) const
    {
        os << " Large=";
        Evaluator::Write(os);
        return os;
    }
};

template <class Evaluator> class EOrbital_RKBS_IBS 
    : public Orbital_RKBS_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    EOrbital_RKBS_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,-1,0) //fix kappa=-1, l=0
    {};

    virtual size_t GetNumFunctions() const {return Evaluator::size();}
    virtual std::ostream& Write(std::ostream& os) const
    {
        os << " Small=";
        Evaluator::Write(os);
        return os;
    }
};

template <class LEvaluator, class SEvaluator> class EOrbital_RKB_IBS 
    : public virtual Orbital_RKB_IBS<double>
    , public Orbital_RKB_IBS_Imp<double>
    , public  BasisSet1::IrrepBasisSetImp<double>
{
public:
    EOrbital_RKB_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_RKB_IBS_Imp(
                new EOrbital_RKBL_IBS<LEvaluator>(N,rmin,rmax,yl),
                new EOrbital_RKBS_IBS<SEvaluator>(N,rmin,rmax,yl)
            )
    , BasisSet1::IrrepBasisSetImp<double>(yl)
    {};

    virtual size_t GetNumFunctions() const {return Orbital_RKB_IBS_Imp<double>::GetNumFunctions();}
};


// Full basis set.
template <class Evaluator> class BasisSet
    : public virtual ::BasisSet1::BasisSet<double>
    , public BasisSet1::BasisSetImp<double>
    , public Evaluator 
{
    using oibs_t=Orbital_IBS<typename Evaluator::IBS_Evaluator_t>; //Corresponding Orbital IBS type
public:
    BasisSet(size_t N, double remin, double remax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(this,N,remin,remax,ir));  
     
        Evaluator::BuildCache(LMax);
    }
    BasisSet(const rvec_t& es, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(this,es,ir));  
     
        Evaluator::BuildCache(LMax);
    }
private:
    void Insert(oibs_t* oibs)
    {
        BasisSet1::BasisSetImp<double>::Insert(oibs);
        Evaluator::Register(oibs);
    }
};

template <class LEvaluator, class SEvaluator> class BasisSet_RKB
    : public virtual ::BasisSet1::BasisSet<double>
    , public BasisSet1::BasisSetImp<double>
{
    using oibs_t=EOrbital_RKB_IBS<LEvaluator,SEvaluator>;
public:
    BasisSet_RKB(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(N,rmin,rmax,ir));  
     
    }
private:
    void Insert(oibs_t* oibs)
    {
        BasisSet1::BasisSetImp<double>::Insert(oibs);
        // Evaluator::Register(oibs->GetEvaluator());
    }
};


}} //namespaces 


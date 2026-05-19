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
import qchem.BasisSet1.DB_Cache;

export 
namespace BasisSet1 {
namespace Atom {



template <isFull_NR_Evaluator Evaluator> class Orbital_IBS 
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

};

template <isFull_NR_Evaluator Evaluator> class Orbital_HF2_IBS 
    : public Orbital_1E_IBS<Evaluator>
    , public Orbital_DFT_IBS<Evaluator>
    , public Orbital_HF1_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    Orbital_HF2_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,yl)
    {};
    Orbital_HF2_IBS(const rvec_t& es, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
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

};

template <is1E_Evaluator Evaluator> class Orbital_1E_HF_IBS 
    : public Orbital_1E_IBS<Evaluator>
    , public Orbital_HF_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    Orbital_1E_HF_IBS(BS_Evaluator* bse,size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_HF_IBS<Evaluator>(bse)
    , IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,yl)
    {};
    Orbital_1E_HF_IBS(BS_Evaluator* bse,const rvec_t& es, const Irrep_QNs::sym_t& yl)
    : Orbital_HF_IBS<Evaluator>(bse)
    , IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(es,yl)
    {};

};

template <is1E_HF_Evaluator Evaluator> class Orbital_1E_HF2_IBS 
    : public Orbital_1E_IBS<Evaluator>
    , public Orbital_HF1_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    Orbital_1E_HF2_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,yl)
    {};
    Orbital_1E_HF2_IBS(const rvec_t& es, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(es,yl)
    {};


   

};

template <isRKBL_Evaluator Evaluator> class EOrbital_RKBL_IBS 
    : public Orbital_RKBL_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    EOrbital_RKBL_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,yl)
    {};


    virtual std::ostream& Write(std::ostream& os) const
    {
        os << " Large=";
        Evaluator::Write(os);
        return os;
    }
};

template <is1E_Evaluator Evaluator> class EOrbital_RKBS_IBS 
    : public Orbital_RKBS_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
public:
    EOrbital_RKBS_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : IrrepBasisSetImp<Evaluator>(yl)
    , Evaluator(N,rmin,rmax,-1,0) //fix kappa=-1, l=0
    {};

    virtual std::ostream& Write(std::ostream& os) const
    {
        os << " Small=";
        Evaluator::Write(os);
        return os;
    }
};

template <isRKBL_Evaluator LEvaluator, is1E_Evaluator SEvaluator> class EOrbital_RKB_IBS 
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

template <class Evaluator> class BasisSet_1E_HF
    : public virtual ::BasisSet1::BasisSet<double>
    , public BasisSet1::BasisSetImp<double>
    , public Evaluator 
{
    using oibs_t=Orbital_1E_HF_IBS<typename Evaluator::IBS_Evaluator_t>; //Corresponding Orbital IBS type
public:
    BasisSet_1E_HF(size_t N, double remin, double remax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(this,N,remin,remax,ir));  
     
        Evaluator::BuildCache(LMax);
    }
    BasisSet_1E_HF(const rvec_t& es, const ElectronConfiguration& ec)
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

template <class Evaluator> class BasisSet_HF2
    : public virtual ::BasisSet1::BasisSet<double>
    , public BasisSet1::BasisSetImp<double>
{
    using oibs_t=Orbital_HF2_IBS<Evaluator>; //Corresponding Orbital IBS type
public:
    BasisSet_HF2(size_t N, double remin, double remax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(N,remin,remax,ir));  
     
    }
    BasisSet_HF2(const rvec_t& es, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(es,ir));  
     
    }
private:
    void Insert(oibs_t* oibs)
    {
        BasisSet1::theGlobalCache->Register(oibs);
        BasisSet1::BasisSetImp<double>::Insert(oibs);
    }
};

template <class Evaluator> class BasisSet_1E_HF2
    : public virtual ::BasisSet1::BasisSet<double>
    , public BasisSet1::BasisSetImp<double>
{
    using oibs_t=Orbital_1E_HF2_IBS<Evaluator>; //Corresponding Orbital IBS type
public:
    BasisSet_1E_HF2(size_t N, double remin, double remax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(N,remin,remax,ir));  
     
    }
    BasisSet_1E_HF2(const rvec_t& es, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(es,ir));  
     
    }
private:
    void Insert(oibs_t* oibs)
    {
        BasisSet1::theGlobalCache->Register(oibs);
        BasisSet1::BasisSetImp<double>::Insert(oibs);
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
        // BasisSet1::theGlobalCache->Register(oibs);
        BasisSet1::BasisSetImp<double>::Insert(oibs);
        // Evaluator::Register(oibs->GetEvaluator());
    }
};


}} //namespaces 


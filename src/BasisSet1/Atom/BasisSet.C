// File: BasisSet/Atom/BSpline/NR/BSpline_BS_Evaluator.C BSpline Basis Set for atoms.
module;
#include <blaze/Math.h>
export module qchem.BasisSet.Atom.BasisSet;

export import qchem.BasisSet;
export import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Fit_IBS;

export import qchem.Symmetry.AtomEC;
export import qchem.Symmetry.Irrep;

import qchem.Symmetry.Yl;
import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.DB_Cache;

export 
namespace BasisSet {
namespace Atom {




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


    virtual ::BasisSet::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const 
    {
        return new Fit_IBS(Evaluator::Rescale(2.0));
    }
    virtual ::BasisSet::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const
    {
        return new Fit_IBS(Evaluator::Rescale(2.0/3.0));
    }

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
    , public ::BasisSet::IrrepBasisSetImp<double>
{
public:
    EOrbital_RKB_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_RKB_IBS_Imp(
                new EOrbital_RKBL_IBS<LEvaluator>(N,rmin,rmax,yl),
                new EOrbital_RKBS_IBS<SEvaluator>(N,rmin,rmax,yl)
            )
    , ::BasisSet::IrrepBasisSetImp<double>(yl)
    {};

    virtual size_t GetNumFunctions() const {return Orbital_RKB_IBS_Imp<double>::GetNumFunctions();}
};


// Full basis set.


template <class Evaluator> class BasisSet_HF2
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
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
        ::BasisSet::theGlobalCache->Register(oibs);
        ::BasisSet::BasisSetImp<double>::Insert(oibs);
    }
};

template <class Evaluator> class BasisSet_1E_HF2
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
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
        ::BasisSet::theGlobalCache->Register(oibs);
        ::BasisSet::BasisSetImp<double>::Insert(oibs);
    }
};

template <class LEvaluator, class SEvaluator> class BasisSet_RKB
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
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
        ::BasisSet::BasisSetImp<double>::Insert(oibs);
    }
};


}} //namespaces 


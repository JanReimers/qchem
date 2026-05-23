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
// import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.BasisSet.Atom.Evaluators.Concepts;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.Internal.DB_Cache;

export 
namespace BasisSet {
namespace Atom {

using namespace Evaluators;

template <isFull_NR_Evaluator Evaluator> class BasisSet_HF
    : public virtual Real_BS
    , public BasisSetImp<double>
{
public:
    // Define the Orbital IBS using mixins.
    class EOrbital_HF_IBS 
        : public Orbital_1E_IBS<Evaluator>
        , public Orbital_DFT_IBS<Evaluator>
        , public Orbital_HF_IBS<Evaluator>
        , public IrrepBasisSetImp<Evaluator>
        , public Evaluator
    {
    public:
        EOrbital_HF_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
        : IrrepBasisSetImp<Evaluator>(yl)
        , Evaluator(N,rmin,rmax,yl)
        {
            theGlobalCache->Register(this); //Can this move to the evaluator level?
        };
        EOrbital_HF_IBS(const rvec_t& es, const Irrep_QNs::sym_t& yl, size_t ltrim=0)
        : IrrepBasisSetImp<Evaluator>(yl)
        , Evaluator(es,yl,ltrim)
        {
            theGlobalCache->Register(this);
        };


        virtual Fit_IBS* CreateCDFitBasisSet(const Cluster*) const 
        {
            return new EFit_IBS(Evaluator::Rescale(2.0));
        }
        virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const
        {
            return new EFit_IBS(Evaluator::Rescale(2.0/3.0));
        }

    };

    BasisSet_HF(size_t N, double remin, double remax, const Atom_EC& ec)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new EOrbital_HF_IBS(N,remin,remax,ir));  
     
    }
    BasisSet_HF(const rvec_t& es, const Atom_EC& ec, size_t ltrim=0)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new EOrbital_HF_IBS(es,ir,ltrim));  
     
    }
};

template <is1E_HF_Evaluator Evaluator> class BasisSet_1E_HF
    : public virtual Real_BS
    , public BasisSetImp<double>
{
public:
    // Define the Orbital IBS using mixins.
    class Orbital_1E_HF_IBS 
    : public Orbital_1E_IBS<Evaluator>
    , public Orbital_HF_IBS<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
    {
    public:
        Orbital_1E_HF_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
        : IrrepBasisSetImp<Evaluator>(yl)
        , Evaluator(N,rmin,rmax,yl)
        {
            theGlobalCache->Register(this);
        };
        Orbital_1E_HF_IBS(const rvec_t& es, const Irrep_QNs::sym_t& yl)
        : IrrepBasisSetImp<Evaluator>(yl)
        , Evaluator(es,yl)
        {
            theGlobalCache->Register(this);
        };
    };

    BasisSet_1E_HF(size_t N, double remin, double remax, const Atom_EC& ec)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new Orbital_1E_HF_IBS(N,remin,remax,ir));  
     
    }
    BasisSet_1E_HF(const rvec_t& es, const Atom_EC& ec)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new Orbital_1E_HF_IBS(es,ir));  
     
    }

};

template <class LEvaluator, class SEvaluator> class BasisSet_RKB
    : public virtual Real_BS
    , public BasisSetImp<double>
{
public:
    class EOrbital_RKBL_IBS 
        : public Orbital_RKBL_IBS<LEvaluator>
        , public IrrepBasisSetImp<LEvaluator>
        , public LEvaluator
    {
    public:
        EOrbital_RKBL_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
        : IrrepBasisSetImp<LEvaluator>(yl)
        , LEvaluator(N,rmin,rmax,yl)
        {};


        virtual std::ostream& Write(std::ostream& os) const
        {
            os << " Large=";
            LEvaluator::Write(os);
            return os;
        }
    };
    class EOrbital_RKBS_IBS 
        : public Orbital_RKBS_IBS<SEvaluator>
        , public IrrepBasisSetImp<SEvaluator>
        , public SEvaluator
    {
    public:
        EOrbital_RKBS_IBS(size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
        : IrrepBasisSetImp<SEvaluator>(yl)
        , SEvaluator(N,rmin,rmax,-1,0) //fix kappa=-1, l=0
        {};

        virtual std::ostream& Write(std::ostream& os) const
        {
            os << " Small=";
            SEvaluator::Write(os);
            return os;
        }
    };
    class EOrbital_RKB_IBS 
        : public virtual Orbital_RKB_IBS<double>
        , public Orbital_RKB_IBS_Imp<double>
        , public ::BasisSet::IrrepBasisSetImp<double>
    {
    public:
        EOrbital_RKB_IBS(size_t N, double remin, double remax, const Irrep_QNs::sym_t& yl)
        : Orbital_RKB_IBS_Imp(
                    new EOrbital_RKBL_IBS(N,remin,remax,yl),
                    new EOrbital_RKBS_IBS(N,remin,remax,yl)
                )
        , IrrepBasisSetImp<double>(yl)
        {};

        virtual size_t GetNumFunctions() const {return Orbital_RKB_IBS_Imp<double>::GetNumFunctions();}
    };
    BasisSet_RKB(size_t N, double ermin, double ermax, const Atom_EC& ec)
    {
        for (auto ir:ec.GetIrreps())
            Insert(new EOrbital_RKB_IBS(N,ermin,ermax,ir));  
     
    }

};


}} //namespaces 


// File: BasisSet/Atom/BSpline/NR/BSpline_BS.C BSpline Basis Set for atoms.
module;
export module qchem.BasisSet1.Atom.BasisSet;

export import qchem.BasisSet1;
export import qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.Symmetry.AtomEC;
export import qchem.Symmetry.Irrep;

import qchem.BasisSet.Atom.BS_Evaluator;
import qchem.BasisSet.Atom.IBS_Evaluator;
import qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Internal.Common;

export 
namespace BasisSet1 {
namespace Atom {

template <class Evaluator> class Orbital_IBS
    : public Orbital_HF_IBS
    , private Evaluator
{
public:
    Orbital_IBS(BS_Evaluator* bse,size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_HF_IBS(bse,yl)
    , Evaluator(N,rmin,rmax,yl)
    {};

    // virtual ::Fit_IBS* CreateCDFitBasisSet(const Real_BS*,const Cluster*) const {return 0;}
    // virtual ::Fit_IBS* CreateVxcFitBasisSet(const Real_BS*,const Cluster*) const {return 0;}
    virtual size_t GetNumFunctions() const {return Evaluator::size();}
    virtual const IBS_Evaluator* GetEvaluator() const {return this;}
    virtual       IBS_Evaluator* GetEvaluator()       {return this;}
};

// template <size_t K> class Fit_IBS 
// : public BSpline_IBS<K>
// , public AtomBS::IrrepBasisSet1
// // , public AtomBS::Fit_IBS
// {
// public:
//     Fit_IBS(size_t N, double rmin, double rmax, size_t L);
// };

    // Full basis set.

template <class Evaluator> class BasisSet
    : public virtual ::BasisSet1::BasisSet<double>
    , public BasisSet1::BS_Common<double>
    , public Evaluator 
{
    using oibs_t=Orbital_IBS<typename Evaluator::IBS_Evaluator_t>; //Corresponding Orbital IBS type
public:
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(this,N,rmin,rmax,ir));  
     
        Evaluator::BuildCache(LMax);
    }
private:
    void Insert(oibs_t* oibs)
    {
        BasisSet1::BS_Common<double>::Insert(oibs);
        Evaluator::Register(oibs->GetEvaluator());
    }
};


}} //namespaces 


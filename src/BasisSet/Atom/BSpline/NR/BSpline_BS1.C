// File: BasisSet/Atom/BSpline/NR/BSpline_BS.C BSpline Basis Set for atoms.
module;
export module qchem.BasisSet.Atom.BSpline.NR.BS1;

export import qchem.BasisSet1;
import qchem.BasisSet.Atom.BS_Evaluator;
import qchem.BasisSet.Internal.Common1;

export import qchem.Orbital_HF_IBS1;
import qchem.BasisSet.Atom.IBS1;

export import qchem.Symmetry.AtomEC;
export import qchem.Symmetry.Irrep;

export namespace AtomBS
{
namespace BSpline1
{
template <size_t K,class Evaluator> class Orbital_IBS
    : public AtomBS::Orbital_HF_IBS1
    , private Evaluator
{
public:
    Orbital_IBS(BS_Evaluator* bse,size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : AtomBS::Orbital_HF_IBS1(bse,yl)
    , Evaluator(N,rmin,rmax,yl)
    {};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet1*,const Cluster*) const {return 0;}
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet1*,const Cluster*) const {return 0;}
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

template <size_t K,template<size_t> class Evaluator> class BasisSet
    : public virtual ::BasisSet1
    , public ::BS_Common1
    , public Evaluator<K> 
{
    using oibs_t=Orbital_IBS<K,typename Evaluator<K>::IBS_Evaluator_t>; //Corresponding Orbital IBS type
public:
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
    {
        const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
        size_t LMax=aec.GetLMax();
        for (auto ir:aec.GetIrreps())
            Insert(new oibs_t(this,N,rmin,rmax,ir));  
     
        Evaluator<K>::BuildCache(LMax);
    }
private:
    void Insert(oibs_t* oibs)
    {
        ::BS_Common1::Insert(oibs);
        Evaluator<K>::Register(oibs->GetEvaluator());
    }
};


}} //namespace AtomBS::BSpline1


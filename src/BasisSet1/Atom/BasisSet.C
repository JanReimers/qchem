// File: BasisSet/Atom/BSpline/NR/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <blaze/Math.h>
export module qchem.BasisSet1.Atom.BasisSet;

export import qchem.BasisSet1;
export import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet1.Fit_IBS;
import qchem.BasisSet1.Atom.IE;

export import qchem.Symmetry.AtomEC;
export import qchem.Symmetry.Irrep;

import qchem.Symmetry.Yl;
import qchem.BasisSet.Atom.BS_Evaluator;
import qchem.BasisSet.Atom.IBS_Evaluator;
import qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Internal.Common;

export 
namespace BasisSet1 {
namespace Atom {

template <class Evaluator> class Fit_IBS
    : public virtual BasisSet1::Fit_IBS 
    , public virtual Integrals_Base
    , public IrrepBasisSetImp
    , public Evaluator
{
    public:
    Fit_IBS(const Evaluator& e) : IrrepBasisSetImp(Irrep_QNs::sym_t(new Yl_Sym(0))), Evaluator(e) {};

    virtual size_t GetNumFunctions() const {return Evaluator::size();}
    virtual const IBS_Evaluator* GetEvaluator() const {return this;}
    virtual       IBS_Evaluator* GetEvaluator()       {return this;}
    
    virtual rsmat_t MakeOverlap  (                ) const {return GetEvaluator()->Overlap   ( );}
    virtual rsmat_t MakeRepulsion(                ) const {return GetEvaluator()->Repulsion ( );}
    virtual  rmat_t MakeRepulsion(const BasisSet1::Fit_IBS& f) const 
    {
        return GetEvaluator()->XRepulsion(dynamic_cast<const IBS_Evaluator&>(f));
    }
    virtual  rvec_t MakeCharge   (                ) const {return GetEvaluator()->Charge    ( );}
    virtual std::ostream&  Write(std::ostream& os) const
    {
        os << "Atom fit IBS ";
        BasisSet1::IrrepBasisSet<double>::Write(os);
        Evaluator::Write(os);
        return os;
    }

};


template <class Evaluator> class Orbital_IBS 
    : public Orbital_1E_IBS
    , public Orbital_DFT_IBS
    , public Orbital_HF_IBS
    , public IrrepBasisSetImp
    , public Evaluator
{
public:
    Orbital_IBS(BS_Evaluator* bse,size_t N, double rmin, double rmax, const Irrep_QNs::sym_t& yl)
    : Orbital_HF_IBS(bse)
    , IrrepBasisSetImp(yl)
    , Evaluator(N,rmin,rmax,yl)
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
    virtual const IBS_Evaluator* GetEvaluator() const {return this;}
    virtual       IBS_Evaluator* GetEvaluator()       {return this;}
};

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


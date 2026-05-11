// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include <cassert>
#include "forward.H"

export module qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Atom.IE;
import qchem.BasisSet1.Internal.IrrepBasisSetImp;
import qchem.BasisSet1.Internal.Orbital_DHF_IBS;
import qchem.BasisSet1.IrrepBasisSet;
import qchem.BasisSet1.Orbital_1E_IBS;
import qchem.BasisSet1.Orbital_DFT_IBS;
import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet1.Atom.BS_Evaluator;

export namespace BasisSet1
{
namespace Atom
{
//
//  Common IrrepBasisSet functionality for atom basis sets.  All the work is done by the evaluator
//
class IrrepBasisSetImp
    : public virtual BasisSet1::IrrepBasisSet<double>
    , public virtual IrrepBasisSet_IDs
    , public virtual Integrals_Base
    , public BasisSet1::IrrepBasisSetImp<double> //Pulls in Symmetry support
{
public:
    IrrepBasisSetImp(const Irrep_QNs::sym_t& yl) : BasisSet1::IrrepBasisSetImp<double>(yl) {}
    virtual size_t  GetNumFunctions() const {return GetEvaluator()->size();}
    // virtual size_t  size           () const {return itsEval->size();}
    virtual rvec_t     operator() (const rvec3_t& r) const {return GetEvaluator()->operator()(r);}
    virtual rvec3vec_t Gradient   (const rvec3_t& r) const {return GetEvaluator()->Gradient(r);}
    // virtual std::ostream&  Write(std::ostream& os) const
    // {
    //     os << GetEvaluator()->Name() << " ";
    //     GetEvaluator()->Write(os);
    //     return os;
    // }
   
    virtual std::string RadialID () const {return GetEvaluator()->RadialID();}
    virtual std::string AngularID() const {return GetEvaluator()->AngularID();}
    virtual std::string Name     () const {return GetEvaluator()->Name();}
};

//
//  1E orbital for atoms.  Use mixins to get the integral evaluations.
//
class Orbital_1E_IBS
    : public virtual BasisSet1::Orbital_1E_IBS<double> //This part has the symmetry.
    , private Integrals_Overlap
    , private Integrals_Kinetic
    , private Integrals_Nuclear
{
public:
    virtual std::ostream&  Write(std::ostream& os) const
    {
        os << "Orbital IBS " << Name() << " ";
        os << "Symmetry=" << GetSymmetry() << " ";
        GetEvaluator()->Write(os);
        return os;
    }
};


class Orbital_DFT_IBS
    : public virtual Integrals_Base
    , public virtual BasisSet1::Orbital_DFT_IBS<double>
{
protected:
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const
    {
        return GetEvaluator()->Overlap(dynamic_cast<const ::IBS_Evaluator&>(c));
    }
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const
    {
        return GetEvaluator()->Repulsion(dynamic_cast<const ::IBS_Evaluator&>(c));
    }
};



class Orbital_HF_IBS
    : public virtual BasisSet1::Orbital_HF_IBS<double> 
    , public virtual Integrals_Base

{
protected:
    Orbital_HF_IBS(BS_Evaluator* bse)  : itsEvaluator(bse) {assert(itsEvaluator);} 

    virtual ERI4 MakeDirect  (const BasisSet1::Orbital_HF_IBS<double>& c) const 
    {
        return itsEvaluator->Direct(GetEvaluator(),dynamic_cast<const Orbital_HF_IBS&>(c).GetEvaluator());
    }
    virtual ERI4 MakeExchange(const BasisSet1::Orbital_HF_IBS<double>& c) const 
    {
        return itsEvaluator->Exchange(GetEvaluator(),dynamic_cast<const Orbital_HF_IBS&>(c).GetEvaluator());
    }
private: 
    BS_Evaluator* itsEvaluator;
};

// Orbital_RKB_IBS does all its integrals by by combining blocks from the L/S sectors.  
// class Orbital_RKB_IBS
//     : public virtual BasisSet1::Orbital_DHF_IBS<double> 
//     , public virtual Integrals_Base
//     // , public Orbital_1E_IBS //pick O
// {
//     virtual       smat_t<T>  MakeRestMass() const;  
// }

class Orbital_RKBL_IBS
    : public virtual BasisSet1::Orbital_RKBL_IBS<double> 
    , public virtual Integrals_Base
    , private Integrals_Overlap
    , private Integrals_Nuclear
{
public:
    virtual rmat_t  MakeKinetic(const Orbital_RKBS_IBS<double>* rkbs) const
    {
        return GetEvaluator()->XKinetic(dynamic_cast<const ::IBS_Evaluator*>(rkbs));
    }
    
};

class Orbital_RKBS_IBS
    : public virtual BasisSet1::Orbital_RKBS_IBS<double> 
    , public virtual Integrals_Base
    , private Integrals_Kinetic
    , private Integrals_Nuclear //RKBS Evaluator overrides Inv_r1 definition
{
    virtual rsmat_t MakeOverlap() const
    {
        return MakeKinetic();
    }
  

};




}} //namespaces
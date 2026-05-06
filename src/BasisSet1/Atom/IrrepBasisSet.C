// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include <cassert>
#include "forward.H"

export module qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Atom.IE;
import qchem.BasisSet1.Internal.IrrepBasisSetImp;
import qchem.BasisSet1.IrrepBasisSet;
import qchem.BasisSet1.Orbital_1E_IBS;
import qchem.BasisSet1.Orbital_DFT_IBS;
import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet.Atom.BS_Evaluator;

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
    , public virtual Streamable
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
//  1E orbital for atoms.  Use mixins to get he integral evaluations.
//
class Orbital_1E_IBS
    : public BasisSet1::Orbital_1E_IBS<double> //This part has the symmetry.
    , private Integrals_Overlap
    , private Integrals_Kinetic
    , private Integrals_Nuclear
{
public:
    virtual std::ostream&  Write(std::ostream& os) const
    {
        os << "Orbital IBS " << Name() << " ";
        BasisSet1::IrrepBasisSet<double>::Write(os);
        GetEvaluator()->Write(os);
        return os;
    }
};


class Orbital_DFT_IBS
    : public virtual Integrals_Base
    , public BasisSet1::Orbital_DFT_IBS<double>
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
    : public BasisSet1::Orbital_HF_IBS<double> 
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

// // Orbital_RKB_IBS does all its integrals in BasisSet.Orbital_RKB_IBS_Common by 
// // by combining blocks from the L/S sectors.  So we just declate RKBL/RKBS here.

// template <class T> class Orbital_RKBL_IBS
//     : public Orbital_RKBL_IBS_Common<T>
//     , public AtomIE_RKBL<T>
// {
// protected:
//     Orbital_RKBL_IBS(const DB_cache<T>* db,const IBS_Evaluator* eval,int kappa)
//         : Orbital_RKBL_IBS_Common<T>(kappa)
//         , AtomIE_RKBL<T>(db,eval)
//         {}
// };

// template <class T> class Orbital_RKBS_IBS
//     : public Orbital_RKBS_IBS_Common<T>
//     , public AtomIE_RKBS<T>
// {
//     virtual const smat_t<T>& Overlap() const {return DB_Kinetic<T>::Kinetic();}
// protected:
//     Orbital_RKBS_IBS(const DB_cache<T>* db,const IBS_Evaluator* eval,int kappa)
//         : Orbital_RKBS_IBS_Common<T>(kappa)
//         , AtomIE_RKBS<T>(db,eval)
//         {}

// };




}} //namespaces
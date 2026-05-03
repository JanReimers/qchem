// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include "forward.H"

export module qchem.BasisSet1.Atom.IBS;
import qchem.BasisSet1.Atom.IE;
import qchem.IrrepBasisSet1;
import qchem.Orbital_1E_IBS1;
import qchem.Orbital_HF_IBS1;
import qchem.BasisSet.Atom.BS_Evaluator;

export namespace AtomBS
{

//
//  Common IrrepBasisSet functionality for atom basis sets.  All the work is done by the evaluator
//
class IrrepBasisSet1
    : public virtual VectorFunction<double>
    , public virtual IrrepBasisSet_IDs
    , public virtual Streamable
    , public virtual Integrals_Base
{
public:
    virtual size_t  GetNumFunctions() const {return GetEvaluator()->size();}
    // virtual size_t  size           () const {return itsEval->size();}
    virtual rvec_t     operator() (const rvec3_t& r) const {return GetEvaluator()->operator()(r);}
    virtual rvec3vec_t Gradient   (const rvec3_t& r) const {return GetEvaluator()->Gradient(r);}
    virtual std::ostream&  Write(std::ostream& os) const
    {
        os << GetEvaluator()->Name() << " ";
        GetEvaluator()->Write(os);
        return os;
    }

    virtual std::string RadialID () const {return GetEvaluator()->RadialID();}
    virtual std::string AngularID() const {return GetEvaluator()->AngularID();}
    virtual std::string Name     () const {return GetEvaluator()->Name();}
};

//
//  1E orbital for atoms.  Use mixins to get he integral evaluations.
//
class Orbital_1E_IBS1
    : public ::Orbital_1E_IBS1<double> //This part has the symmetry.
    , public AtomBS::IrrepBasisSet1
    , private AtomBS::Integrals_Overlap1
    , private AtomBS::Integrals_Kinetic1
    , private AtomBS::Integrals_Nuclear1
{
public:
    Orbital_1E_IBS1(const Irrep_QNs::sym_t& yl) : ::Orbital_1E_IBS1<double>(yl) {};

};


// template <class T> class Orbital_DFT_IBS
//     : public virtual ::Orbital_DFT_IBS<T> 
//     , public Orbital_DFT_IBS_Common<double>
//     , public AtomIE_DFT<double>
// {
// protected:
//     Orbital_DFT_IBS(const DB_cache<double>* db,const IBS_Evaluator* eval) 
//     : AtomIE_DFT<double>(db,eval)
//     {};
// };

class Orbital_HF_IBS1
    : public ::Orbital_HF_IBS1<double> 
    , public AtomBS::Orbital_1E_IBS1
    , public AtomBS::Integrals_HF1
{
protected:
    Orbital_HF_IBS1(BS_Evaluator* bse,const Irrep_QNs::sym_t& yl) 
        : AtomBS::Orbital_1E_IBS1(yl) 
        , AtomBS::Integrals_HF1(bse)
        {} 
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

// class Fit_IBS
//     : public virtual ::Fit_IBS 
//     , public Fit_IBS_Common
//     , public AtomIE_Overlap<double>
//     , public AtomIE_Fit
// {
// protected:
//     Fit_IBS(const DB_cache<double>* db,const IBS_Evaluator* eval) 
//     : AtomIE_Overlap<double>(db,eval)
//     , AtomIE_Fit(db,eval)
//     {};
// };



} //namespace
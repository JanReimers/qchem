// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include "forward.H"

export module qchem.BasisSet.Atom.IBS1;
import qchem.BasisSet.Atom.IE1;
import qchem.IrrepBasisSet1;
import qchem.BasisSet.DB_Cache1;

export namespace AtomBS
{

class IrrepBasisSet1
    : public virtual VectorFunction<double>
    , public virtual IrrepBasisSet_IDs
    , public virtual Streamable
{
public:
    IrrepBasisSet1(IBS_Evaluator* eval)
    : itsEval(eval)
    {};
    virtual size_t  GetNumFunctions() const {return itsEval->size();}
    // virtual size_t  size           () const {return itsEval->size();}
    virtual rvec_t     operator() (const rvec3_t& r) const {return itsEval->operator()(r);}
    virtual rvec3vec_t Gradient   (const rvec3_t& r) const {return itsEval->Gradient(r);}
    virtual std::ostream&  Write(std::ostream& os) const
    {
        os << itsEval->Name() << " ";
        itsEval->Write(os);
        return os;
    }

    virtual std::string RadialID () const {return itsEval->RadialID();}
    virtual std::string AngularID() const {return itsEval->AngularID();}
    virtual std::string Name     () const {return itsEval->Name();}
private:
    IBS_Evaluator* itsEval;
};

// template <class T> class Orbital_IBS1
//     : public virtual ::Orbital_IBS1<T> 
//     , public AtomIE_Overlap1<double>
//     , public AtomIE_Kinetic1<double>
//     , public AtomIE_Nuclear1<double>
// {
// protected:
//     Orbital_IBS1(const IBS_Evaluator* eval) 
//     : AtomIE_Overlap1<double>(eval)
//     , AtomIE_Kinetic1<double>(eval)
//     , AtomIE_Nuclear1<double>(eval) 
//     {};
// };

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

// template <class T> class Orbital_HF_IBS
//     : public virtual ::Orbital_HF_IBS<T> 
//     , public Orbital_HF_IBS_Common<T>
// {
// protected:
//     Orbital_HF_IBS(const Integrals_BS_HF<double>* db) : Orbital_HF_IBS_Common<T>(db) {};
// };

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
// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
export module qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.IrrepBasisSet;

export namespace Atom
{

class IrrepBasisSet
    : public virtual Real_IBS
    , public         IrrepBasisSet_Common<double>
{
    typedef typename VectorFunction<double>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<double>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
public:
    IrrepBasisSet(IBS_Evaluator* eval, Symmetry* sym)
    : IrrepBasisSet_Common<double>(sym)
    , itsEval(eval)
    {};
    virtual size_t  GetNumFunctions() const {return itsEval->size();}
    virtual size_t  size           () const {return itsEval->size();}
    virtual Vec     operator() (const RVec3& r) const {return itsEval->operator()(r);}
    virtual Vec3Vec Gradient   (const RVec3& r) const {return itsEval->Gradient(r);}
    virtual std::ostream&  Write(std::ostream&) const;
private:
    IBS_Evaluator* itsEval;
};

template <class T> class Orbital_IBS
    : public virtual ::Orbital_IBS<T> //brings in Integrals_Overlap<T>
    , public Orbital_IBS_Common<double>
    , public AtomIE_Overlap<double>
    , public AtomIE_Kinetic<double>
    , public AtomIE_Nuclear<double>
{
protected:
    Orbital_IBS(const DB_cache<double>* db,const IBS_Evaluator* eval) 
    : AtomIE_Overlap<double>(db,eval)
    , AtomIE_Kinetic<double>(db,eval)
    , AtomIE_Nuclear<double>(db,eval) 
    {};
};

template <class T> class Orbital_DFT_IBS
    : public virtual ::Orbital_DFT_IBS<T> 
    , public Orbital_DFT_IBS_Common<double>
    , public AtomIE_DFT<double>
{
protected:
    Orbital_DFT_IBS(const DB_cache<double>* db,const IBS_Evaluator* eval) 
    : AtomIE_DFT<double>(db,eval)
    {};
};

template <class T> class Orbital_HF_IBS
    : public virtual ::Orbital_HF_IBS<T> 
    , public Orbital_HF_IBS_Common<T>
{
protected:
    Orbital_HF_IBS(const DB_BS_2E<double>* db) : Orbital_HF_IBS_Common<T>(db) {};
};

// Orbital_RKB_IBS does all its integrals in BasisSet.Orbital_RKB_IBS_Common by 
// by combining blocks from the L/S sectors.  So we just declate RKBL/RKBS here.

template <class T> class Orbital_RKBL_IBS
    : public Orbital_RKBL_IBS_Common<T>
    , public AtomIE_RKBL<T>
{
protected:
    Orbital_RKBL_IBS(const DB_cache<T>* db,const IBS_Evaluator* eval,int kappa)
        : Orbital_RKBL_IBS_Common<T>(kappa)
        , AtomIE_RKBL<T>(db,eval)
        {}
};

template <class T> class Orbital_RKBS_IBS
    : public Orbital_RKBS_IBS_Common<T>
    , public AtomIE_RKBS<T>
{
    virtual const SMatrix<T>& Overlap() const {return DB_Kinetic<T>::Kinetic();}
protected:
    Orbital_RKBS_IBS(const DB_cache<T>* db,const IBS_Evaluator* eval,int kappa)
        : Orbital_RKBS_IBS_Common<T>(kappa)
        , AtomIE_RKBS<T>(db,eval)
        {}
};

class Fit_IBS
    : public virtual ::Fit_IBS 
    , public Fit_IBS_Common
    , public AtomIE_Overlap<double>
    , public AtomIE_Fit
{
protected:
    Fit_IBS(const DB_cache<double>* db,const IBS_Evaluator* eval) 
    : AtomIE_Overlap<double>(db,eval)
    , AtomIE_Fit(db,eval)
    {};
};



} //namespace
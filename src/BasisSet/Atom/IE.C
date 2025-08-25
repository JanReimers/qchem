// File: AtomIE.C Common IE code for all atom basis sets.
module;
#include <memory>
#include <cassert>
export module qchem.BasisSet.Atom.IE;
export import qchem.BasisSet.Internal.HeapDB;
export import qchem.BasisSet.Internal.Cache4;
export import oml.Vector;
export import qchem.Orbital_DHF_IBS;
export import qchem.Types;
export import BasisSet.Atom.IBS_Evaluator;
export import BasisSet.Atom.BS_Evaluator;
import qchem.BasisSet.Atom.Internal.BFGrouper;
import qchem.BasisSet.Atom.IEClient;

import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IEClient;
import qchem.Fit_IBS;
import qchem.Orbital_DFT_IBS;

export
{

template <class T> class AtomIE_Overlap
: public DB_Overlap<T>
{
protected:
<<<<<<< HEAD
    AtomIE_Overlap(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_Overlap<T>(db), eval(_eval) {};
=======
    AtomIE_Overlap(const DB_cache<T>* db,const IE_Primatives* _pie,const IBS_Evaluator* _eval) : DB_Overlap<T>(db), pie(_pie), eval(_eval) {};
>>>>>>> origin/main
    virtual SMatrix<T> MakeOverlap() const {return eval->Overlap();}
private:
    const IBS_Evaluator* eval;
};
template <class T> class AtomIE_Kinetic
: public DB_Kinetic<T>
{
protected:
<<<<<<< HEAD
    AtomIE_Kinetic(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_Kinetic<T>(db), eval(_eval) {};
    virtual SMatrix<T> MakeKinetic() const
    {
        int l=eval->Getl();
        return eval->Grad2() + l*(l+1)*eval->Inv_r2();
=======
    AtomIE_Kinetic(const DB_cache<T>* db,const IE_Primatives* _pie,const IBS_Evaluator* _eval) : DB_Kinetic<T>(db), pie(_pie), eval(_eval) {};
    virtual SMatrix<T> MakeKinetic() const
    {
        int l=eval->Getl();
        return eval->Grad2()+l*(l+1)*eval->Inv_r1();
>>>>>>> origin/main
    }
private:
    const IBS_Evaluator* eval;
};
template <class T> class AtomIE_Nuclear
: public DB_Nuclear<T>
{
protected:
    AtomIE_Nuclear(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_Nuclear<T>(db), eval(_eval) {};
    virtual SMatrix<T> MakeNuclear(const Cluster* cl) const
    {
        assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
        return Z*eval->Inv_r1();
    }
private:
    const IBS_Evaluator* eval;
};
template <class T> class AtomIE_XKinetic
: public DB_XKinetic<T>
{
protected:
    virtual Matrix<T> MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const
    {
        return eval->XKinetic(rkbs);
    }
    AtomIE_XKinetic(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_XKinetic<T>(db), eval(_eval) {};
private:
    const IBS_Evaluator* eval;
};

// HF
class AtomIE_BS_2E_Angular
{
public:
    virtual ~AtomIE_BS_2E_Angular() {};
    typedef AtomIrrepIEClient iec_t;
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const=0;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* c) const=0;
};

template <class T> class AtomIE_BS_2E 
    : public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper
{
public:
    AtomIE_BS_2E(AtomIE_BS_2E_Angular* a) : itsAngular(a), itsEvaluator(0) {};
    AtomIE_BS_2E(BS_Evaluator* bse,AtomIE_BS_2E_Angular* a) : itsAngular(a), itsEvaluator(bse) {};
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const IrrepIEClient*);
    virtual void Append(const IrrepIEClient*, IBS_Evaluator*);
private: 
    std::unique_ptr<AtomIE_BS_2E_Angular> itsAngular;
    BS_Evaluator* itsEvaluator;
};

// DFT
template <class T> class AtomIE_DFT 
: public DB_DFT<T>
{
protected:
    AtomIE_DFT(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_DFT<T>(db), eval(_eval) {};
    
    virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const {return eval->Overlap  (c);}
    virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const {return eval->Repulsion(c);}
private:
    const IBS_Evaluator* eval;
};
// DHF
template <class T> class AtomIE_RKBL 
    : public AtomIE_Overlap<T>
    , public AtomIE_XKinetic<T>
    , public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBL(const DB_cache<T>* db,const IBS_Evaluator* eval) 
    : AtomIE_Overlap <T>(db,eval)
    , AtomIE_XKinetic<T>(db,eval)
    , AtomIE_Nuclear <T>(db,eval) 
    {};

};
template <class T> class AtomIE_RKBS 
: public AtomIE_Kinetic<T>
, public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBS(const DB_cache<T>* db,const IBS_Evaluator* eval) 
    : AtomIE_Kinetic<T>(db,eval)
    , AtomIE_Nuclear<T>(db,eval) {};
};
// Fit
class AtomIE_Fit 
: public DB_Fit
{
    protected:
    AtomIE_Fit(const DB_cache<double>* db,const IBS_Evaluator* _eval) : DB_Fit(db), eval(_eval) {};

    virtual  Vector<double> MakeCharge   (                ) const {return eval->Charge    ( );}
    virtual SMatrix<double> MakeRepulsion(                ) const {return eval->Repulsion ( );}
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS& f) const {return eval->XRepulsion(f);}
private:
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
private:
    const IBS_Evaluator* eval;
};

} // export block
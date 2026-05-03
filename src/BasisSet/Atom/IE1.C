// File: BasisSet/Atom/IE.C Common Integral Engine (IE) code for all atom basis sets.
module;
#include <cassert>
#include "blaze/Math.h"
export module qchem.BasisSet.Atom.IE1;
export import qchem.BasisSet.Internal.ERI4;

export import qchem.Orbital_1E_IBS1;
// export import qchem.Orbital_DHF_IBS;
// export import qchem.Orbital_DFT_IBS;
export import qchem.Orbital_HF_IBS1;
// export import qchem.Fit_IBS;

export import qchem.BasisSet.Atom.IBS_Evaluator;
import qchem.BasisSet.Atom.BS_Evaluator;
import qchem.Types;


export namespace AtomBS
{

class Integrals_Base
{
public:
    virtual const IBS_Evaluator* GetEvaluator() const=0;
    virtual IBS_Evaluator* GetEvaluator()=0;
};

class Integrals_Overlap1
: public virtual ::Integrals_Overlap1<double>
, public virtual Integrals_Base
{
protected:
    virtual smat_t<double> MakeOverlap() const {return GetEvaluator()->Overlap();}
};
class Integrals_Kinetic1
: public virtual ::Integrals_Kinetic1<double>
, public virtual Integrals_Base
{
public:
    virtual smat_t<double> MakeKinetic() const
    {
        auto eval=GetEvaluator();
        int l=eval->Getl();
        return eval->Grad2() + l*(l+1)*eval->Inv_r2();
    }
};
class Integrals_Nuclear1
: public virtual ::Integrals_Nuclear1<double>
, public virtual Integrals_Base
{
protected:
    virtual smat_t<double> MakeNuclear(const Cluster* cl) const
    {
        assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
        return Z*GetEvaluator()->Inv_r1();
    }
};

class Integrals_HF1
: public virtual ::Integrals_HF1<double>
, public virtual Integrals_Base
{
protected:
    Integrals_HF1(BS_Evaluator* bse) : itsEvaluator(bse) {assert(itsEvaluator);}
    virtual ERI4 MakeDirect  (const Orbital_HF_IBS1<double>& c) const 
    {
        return itsEvaluator->Direct(GetEvaluator(),dynamic_cast<const Integrals_HF1&>(c).GetEvaluator());
    }
    virtual ERI4 MakeExchange(const Orbital_HF_IBS1<double>& c) const 
    {
        return itsEvaluator->Exchange(GetEvaluator(),dynamic_cast<const Integrals_HF1&>(c).GetEvaluator());
    }
private: 
    BS_Evaluator* itsEvaluator;
};
// template <class T> class AtomIE_XKinetic
// : public DB_XKinetic<T>
// {
// protected:
//     virtual mat_t<T> MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const
//     {
//         return eval->XKinetic(rkbs);
//     }
//     AtomIE_XKinetic(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_XKinetic<T>(db), eval(_eval) {};
// private:
//     const IBS_Evaluator* eval;
// };

// HF

// class AtomIE_BS_HF1 
//     : public Integrals_BS_HF1<double>
// {
// public:
//     AtomIE_BS_HF1(BS_Evaluator* bse) : itsEvaluator(bse) {};
//     using IBS_t=Orbital_HF_IBS1<double>;
//     virtual ERI4 MakeDirect  (const IBS_t* a, const IBS_t* c) const;
//     virtual ERI4 MakeExchange(const IBS_t* a, const IBS_t* c) const;
// protected:
//     virtual void Append(const Orbital_HF_IBS1<T>*, IBS_Evaluator*);

// private: 
//     BS_Evaluator* itsEvaluator;
// };

// DFT
// template <class T> class AtomIE_DFT 
// : public DB_DFT<T>
// {
// protected:
//     AtomIE_DFT(const DB_cache<T>* db,const IBS_Evaluator* _eval) : DB_DFT<T>(db), eval(_eval) {};
    
//     virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const {return eval->Overlap  (c);}
//     virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const {return eval->Repulsion(c);}
// private:
//     const IBS_Evaluator* eval;
// };
// DHF
// template <class T> class AtomIE_RKBL 
//     : public AtomIE_Overlap<T>
//     , public AtomIE_XKinetic<T>
//     , public AtomIE_Nuclear<T>
// {
// protected:
//     AtomIE_RKBL(const DB_cache<T>* db,const IBS_Evaluator* eval) 
//     : AtomIE_Overlap <T>(db,eval)
//     , AtomIE_XKinetic<T>(db,eval)
//     , AtomIE_Nuclear <T>(db,eval) 
//     {};

// };
// template <class T> class AtomIE_RKBS 
// : public AtomIE_Kinetic<T>
// , public AtomIE_Nuclear<T>
// {
// protected:
//     AtomIE_RKBS(const DB_cache<T>* db,const IBS_Evaluator* eval) 
//     : AtomIE_Kinetic<T>(db,eval)
//     , AtomIE_Nuclear<T>(db,eval) {};
// };
// // Fit
// class AtomIE_Fit 
// : public DB_Fit
// {
//     protected:
//     AtomIE_Fit(const DB_cache<double>* db,const IBS_Evaluator* _eval) : DB_Fit(db), eval(_eval) {};

//     virtual  rvec_t MakeCharge   (                ) const {return eval->Charge    ( );}
//     virtual rsmat_t MakeRepulsion(                ) const {return eval->Repulsion ( );}
//     virtual  rmat_t MakeRepulsion(const Fit_IBS& f) const {return eval->XRepulsion(f);}
// private:
//     using DB_Fit::Charge; //un hide
//     using DB_Fit::Repulsion; //un hide
// private:
//     const IBS_Evaluator* eval;
// };

} // export block
// File: BasisSet/Atom/IE.C Common Integral Engine (IE) code for all atom basis sets.
module;
#include <cassert>
#include "blaze/Math.h"
export module qchem.BasisSet1.Atom.IE;
export import qchem.BasisSet1.Internal.ERI4;

export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Atom.IBS_Evaluator;
import qchem.BasisSet1.Atom.BS_Evaluator;
import qchem.Types;


export namespace BasisSet1
{
namespace Atom
{

class Integrals_Base
{
public:
    virtual const IBS_Evaluator* GetEvaluator() const=0;
    virtual IBS_Evaluator* GetEvaluator()=0;
};

//
//  Implement these separately as they shared between NR and RKB 1E orbital IBS implementations.
//
class Integrals_Overlap
: public virtual BasisSet1::Integrals_Overlap<double>
, public virtual Integrals_Base
{
protected:
    virtual smat_t<double> MakeOverlap() const {return GetEvaluator()->Overlap();}
};
class Integrals_Kinetic
: public virtual BasisSet1::Integrals_Kinetic<double>
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
class Integrals_Nuclear
: public virtual BasisSet1::Integrals_Nuclear<double>
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

// class Integrals_HF
// : public virtual BasisSet1::Integrals_HF<double>
// , public virtual Integrals_Base
// {
// protected:
//     Integrals_HF(BS_Evaluator* bse) : itsEvaluator(bse) {assert(itsEvaluator);}
//     virtual ERI4 MakeDirect  (const Orbital_HF_IBS<double>& c) const 
//     {
//         return itsEvaluator->Direct(GetEvaluator(),dynamic_cast<const Integrals_HF&>(c).GetEvaluator());
//     }
//     virtual ERI4 MakeExchange(const Orbital_HF_IBS<double>& c) const 
//     {
//         return itsEvaluator->Exchange(GetEvaluator(),dynamic_cast<const Integrals_HF&>(c).GetEvaluator());
//     }
// private: 
//     BS_Evaluator* itsEvaluator;
// };
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

}} // namespaces
// File: BasisSet/Atom/IE.C Common Integral Engine (IE) code for all atom basis sets.
module;
#include <cassert>
#include "blaze/Math.h"
export module qchem.BasisSet1.Atom.IE;
export import qchem.BasisSet1.Internal.ERI4;

export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Atom.Evaluators.IBS;
import qchem.BasisSet1.Atom.Evaluators.BS;
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
// 

template <isEvaluator E> class Integrals_EOverlap
: public virtual BasisSet1::Integrals_Overlap<double>
, public virtual Integrals_Base
{
protected:
    virtual smat_t<double> MakeOverlap() const 
    {
        auto e=dynamic_cast<const E*>(GetEvaluator());
        size_t N=e->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= e->Overlap(i,j);

        return S;
    }
};

template <isEvaluator E> class Integrals_EKinetic
: public virtual BasisSet1::Integrals_Kinetic<double>
, public virtual Integrals_Base
{
protected:
    virtual smat_t<double> MakeKinetic() const 
    {
        auto e=dynamic_cast<const E*>(GetEvaluator());
        size_t N=e->size();
        int l=e->Getl();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= e->Grad2(i,j) + l*(l+1)*e->Inv_r2(i,j);

        return S;
    }
};

template <isEvaluator E> class Integrals_ENuclear
: public virtual BasisSet1::Integrals_Nuclear<double>
, public virtual Integrals_Base
{
protected:
    virtual smat_t<double> MakeNuclear(const Cluster* cl) const 
    {
        assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
        auto e=dynamic_cast<const E*>(GetEvaluator());
        size_t N=e->size();
        rsmat_t S(N);
        for (auto i:iv_t(0,N))
            for (auto j:iv_t(i,N))
                S(i,j)= Z*e->Inv_r1(i,j);

        return S;
    }
};

}} // namespaces
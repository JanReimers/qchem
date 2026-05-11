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


}} // namespaces
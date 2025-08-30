// File: BasisSet/Atom/radial/Exponential_Evaluator.C  Common for Slater and Gaussian evaluators
module;
export module qchem.BasisSet.Atom.Internal.Exponential_IBS_Evaluator;
export import qchem.BasisSet.Atom.IBS_Evaluator;

export class Exponential_IBS_Evaluator : public IBS_Evaluator
{
public:
    Exponential_IBS_Evaluator(const   ds_t& _es, int l, const is_t& mls) : IBS_Evaluator(l,mls), es(_es ) {};
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.

protected:
    ds_t es; 
};
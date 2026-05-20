// File: BasisSet1/Atom/Evaluators/Internal/Exponential_Evaluator.C  Common base for Slater and Gaussian evaluators
module;
#include <string>
export module qchem.BasisSet.Atom.Evaluators.Internal.Exponential_IBS_Evaluator;
export import qchem.BasisSet.Atom.Evaluators.IBS;

export class Exponential_IBS_Evaluator : public IBS_Evaluator
{
public:
    Exponential_IBS_Evaluator(const   rvec_t& _es, int l, const is_t& mls) 
        : IBS_Evaluator(l,mls), es(_es )
        , isEvenTempered(EvenTempered(es))
        {};
    Exponential_IBS_Evaluator(const   rvec_t& _es, const Irrep_QNs::sym_t& ir) 
        : IBS_Evaluator(ir), es(_es ) 
        , isEvenTempered(EvenTempered(es))
        {};
    virtual ~Exponential_IBS_Evaluator() {}; //g++ 15.2 BUG Compiler implemented destructor not created with -O2.
    virtual void Register(Grouper*); //Set up unique spline or exponent indexes.
    virtual std::string RadialID () const;

protected:
    static bool EvenTempered(const rvec_t&);
    rvec_t es;
    bool isEvenTempered; // e[i]=beta*e[i-1] +/- eps;
};
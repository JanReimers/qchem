// File: BasisSet/Atom/Evaluators/Internal/Exponential_Evaluator.C  Common base for Slater and Gaussian evaluators
module;
#include <string>
#include <vector>
#include <cassert>
#include "forward.H"

export module qchem.BasisSet.Atom.Evaluators.Internal.ExponentialEvaluator;
import qchem.BasisSet.Atom.Evaluators.Internal.Grouper;
export import qchem.BasisSet.Atom.Evaluators;
export import qchem.Symmetry.Spherical;
import qchem.Blaze;

export namespace qchem::BasisSet::Atom::Evaluators
{
    
class ExponentialEvaluator : public virtual Evaluator, public virtual HF_Evaluator
{
public:
    ExponentialEvaluator(const rvec_t& _es, const sym_t& ir, size_t ltrim=0)
        : es()
        , grouper(0)
        , isEvenTempered(EvenTempered(es))
        {
            int l=Symmetry::Getl(ir);
            assert(ltrim<5);
            assert(l<=3);
            size_t nfront=(3-l)*ltrim;
            size_t nback=l*ltrim;
            assert(_es.size()>nfront+nback);
            size_t n=_es.size()-nfront-nback;
            es=blazem::subvector(_es,nfront,n);
        };
    virtual ~ExponentialEvaluator() {}; //g++ 15.2 BUG Compiler implemented destructor not created with -O2.
    virtual void   Register(Grouper*) override; //Set up unique spline or exponent indexes.
    virtual std::string RadialID() const override;
    virtual size_t size() const override { return ns.size(); }
    virtual size_t es_index(size_t i     ) const {return es_indices[i];}

    virtual rvec_t Norm() const override { return ns; }

protected:
    friend ::Cache4Tests;
    static bool EvenTempered(const rvec_t&);
    rvec_t es;
    rvec_t ns;
    const  ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index

    bool isEvenTempered; // e[i]=beta*e[i-1] +/- eps;
};

} //namespace
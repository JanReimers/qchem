// File: BasisSet/Atom/radial/Imp/Evaluator.C
module;
#include <cassert>
#include <blaze/math/DynamicVector.h>

module qchem.BasisSet.Atom.Evaluators.Slater.IBS;
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.Integrals; 
import qchem.Symmetry.Factory;
import qchem.Math;

namespace BasisSet::Atom::Evaluators::Slater
{

RKBS_Evaluator::RKBS_Evaluator(size_t N, double emin, double emax, int _kappa, int l)
    : Evaluator(N,emin,emax,SymmetryFactory::YFactory())
    , kappa(_kappa) 
{
    ns=norms();
}
    
std::string RKBS_Evaluator::Name() const
{
    return "SL RKB ";
}

rvec_t RKBS_Evaluator::norms() const
{
    return Evaluator::norms()/(2*c_light);
}

rvec_t RKBS_Evaluator::eval(const rvec3_t& r) const
{
    double mr=norm(r);
    rvec_t f=-es;
    if (kappa >0) 
        f+=(2*kappa+1)/mr;
        
    return f*slater(mr,l,es,ns);
}

rvec_t RKBS_Evaluator::operator() (const rvec3_t& r) const
{
   return eval(r); 
}

rvec3vec_t RKBS_Evaluator::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    rvec_t grad=eval(r)*(l/mr-es);
   
    // std::cout << "grad=";
    // for (auto g:grad) std::cout << g << " ";
    // std::cout << std::endl;
    rvec3_t rhat=r/norm(r);
    size_t i=0;
    for (auto& g:grad) ret[i++]=g*rhat;
    // std::cout << "ret=";
    // for (auto g:ret) std::cout << g << " ";
    // std::cout << std::endl;
    return ret;
}

} //namespace
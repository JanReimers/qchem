// File: BasisSet/Radial/radial/Imp/Evaluator.C
module;
#include <ostream>   // R1.9: free operator<<(ostream&,const char*) is NOT visible via a module import
#include <cassert>
#include <string>
module qchem.BasisSet.Radial.Evaluators.Slater.IBS;
import qchem.BasisSet.Radial.Evaluators.Slater.Internal.Integrals; 
import qchem.Blaze;

namespace qchem::BasisSet::Radial::Evaluators::Slater
{

std::string RKBS_Evaluator::Name() const
{
    return "SL RKB ";
}

rvec_t RKBS_Evaluator::norms() const
{
    return Radial::norms()/(2*c_light);
}

rvec_t RKBS_Evaluator::eval(const rvec3_t& r) const
{
    int κ = Getκ();
    double mr=norm(r);
    rvec_t f=-es;
    if (κ>0)
        f+=(2*κ+1)/mr;
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
    // d/dr [f(r) g(r)] with f = -ζ (+ (2κ+1)/r for κ>0) and g the large-component Slater: f g' + f' g.  (It was f g' alone until the
    // hand-gradient audit of 2026-11-09.)
    int κ = Getκ();
    rvec_t f=-es, fp(es.size(), 0.0);
    if (κ>0)
    {
        f +=(2*κ+1)/mr;
        fp-=(2*κ+1)/(mr*mr);
    }
    rvec_t grad=f*grad_slater(mr,l,es,ns)+fp*slater(mr,l,es,ns);
   
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
// File: BasisSet/Atom/radial/Imp/Evaluator.C
module;
#include <cassert>
#include <string>
module qchem.BasisSet.Atom.Evaluators.Slater.IBS;
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.Integrals; 
import qchem.Symmetry.Factory;
import qchem.Math;
import qchem.Blaze;

namespace BasisSet::Atom::Evaluators::Slater
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
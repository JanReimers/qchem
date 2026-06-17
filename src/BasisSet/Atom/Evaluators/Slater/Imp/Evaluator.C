// File: BasisSet/Atom/radial/Imp/Evaluator.C
module;
#include <cassert>
#include <iostream>

module qchem.BasisSet.Atom.Evaluators.Slater.IBS;
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.Rk; 
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.Integrals; 
import qchem.BasisSet.Atom.Evaluators.Slater.Internal.ExponentScaler; 
import qchem.Math;
import qchem.Blaze;

namespace BasisSet::Atom::Evaluators::Slater
{

//---------------------------------------------------------------------------
//
//  Start member functions.
//

std::string Radial::Name() const
{
    return "SL";
}
std::string Radial::RadialType() const
{
    return Name();
}

Cache4*    Radial::MakeCache4() const
{
    return new Slater_Cache4();
}


rvec_t Radial::exponents(size_t N, double emin, double emax, const sym_t& ir)
{
    size_t LMax=3; //TODO how do we get the real LMax(Z) into this?
    ::Slater::ExponentScaler ss(N,emin,emax,LMax);
    return ss.Get_es(ir);
}


rvec_t Radial::norms() const
{
    size_t N=es.size();    
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Slater::Integral(2*es[i],2*l)); 
    return ret;
}





rvec_t Radial::operator() (const rvec3_t& r) const
{
    return slater(norm(r),l,es,ns);
}

rvec3vec_t Radial::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    rvec_t grad=grad_slater(norm(r),l,es,ns);
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


std::ostream&  Radial::Write(std::ostream& os) const
{
    return os << " N=" << size() << " α={" << es[0] << " ... " << es[size()-1] << "}";
}

} //namespace
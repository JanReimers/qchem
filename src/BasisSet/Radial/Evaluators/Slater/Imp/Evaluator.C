// File: BasisSet/Radial/radial/Imp/Evaluator.C
module;
#include <cassert>
#include <iostream>
#include <cmath>

module qchem.BasisSet.Radial.Evaluators.Slater.IBS;
import qchem.BasisSet.Radial.Evaluators.Slater.Internal.Rk; 
import qchem.BasisSet.Radial.Evaluators.Slater.Internal.Integrals; 
import qchem.BasisSet.Radial.Evaluators.Slater.Internal.ExponentScaler; 
import qchem.Blaze;
import qchem.IntPow;   // uintpow: a power-m fit function's r^{2m}
import qchem.Math;

namespace qchem::BasisSet::Radial::Evaluators::Slater
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
    ::qchem::Slater::ExponentScaler ss(N,emin,emax,LMax);
    return ss.Get_es(ir);
}


rvec_t Radial::norms() const
{
    size_t N=es.size();    
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::qchem::Slater::Integral(2*es[i], ms.empty() ? 2*l : 4*ms[i]));   // a power-m fit function is r^{2m} e^{-zeta r}
    return ret;
}





rvec_t Radial::operator() (const rvec3_t& r) const
{
    if (ms.empty()) return slater(norm(r),l,es,ns);
    const double mr=norm(r);
    rvec_t v(es.size());
    for (size_t i=0;i<es.size();i++) v[i]=ns[i]*uintpow(mr,2*ms[i])*std::exp(-es[i]*mr);   // a power-m fit function
    return v;
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
    rvec_t grad(size());
    if (ms.empty()) grad=grad_slater(norm(r),l,es,ns);
    else for (size_t i=0;i<es.size();i++) grad[i]=(2*ms[i]/mr-es[i])*ns[i]*uintpow(mr,2*ms[i])*std::exp(-es[i]*mr);   // d/dr r^{2m} e^{-zeta r}
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
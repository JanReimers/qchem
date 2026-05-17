// File: BasisSet/Atom/Gaussian/NR/Imp/IBS_Evaluator.C
module;
#include <cmath>
#include <cassert>
#include <iostream>
#include <blaze/math/DynamicVector.h>

module qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.GaussianIntegrals; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.ExponentScaler; 
import Common.Constants;


//---------------------------------------------------------------------------
//
//  Start member functions.
//
std::string Gaussian_IBS_Evaluator::Name() const
{
    return "Spherical Gaussian ";
}

std::string Gaussian_IBS_Evaluator::RadialType() const
{
    return "SG";
}

Cache41*    Gaussian_IBS_Evaluator::MakeCache4() const
{
    return new Gaussian_Cache4();
}

rvec_t Gaussian_IBS_Evaluator::exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir)
{
    size_t LMax=3; //TODO how do we get the real LMax(Z) into this?
    ::Gaussian::ExponentScaler ss(N,emin,emax,LMax);
    return ss.Get_es(ir);
}

rvec_t Gaussian_IBS_Evaluator::norms() const
{
    size_t N=es.size();    
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(Gaussian::Integral(2*es[i],2*l)); 
    return ret;
}


rvec_t Gaussian_IBS_Evaluator::operator() (const rvec3_t& r) const
{
    return gaussian(norm(r),l,es,ns); 
}

rvec3vec_t Gaussian_IBS_Evaluator::Gradient(const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    double mr=norm(r);
    if (mr==0.0)
    {
        ret=rvec3_t{0,0,0};
        return ret;
    }
    rvec_t grad=grad_gaussian(norm(r),l,es,ns);
    rvec3_t rhat=r/norm(r);
    size_t i=0;
    for (auto& g:grad) ret[i++]=g*rhat;
    return ret;
}

std::ostream&  Gaussian_IBS_Evaluator::Write(std::ostream& os) const
{
    return os << " with N=" << es.size() << " basis functions, alpha={" << es[0] << " ... " << es[size()-1] << "}" << std::endl;
}



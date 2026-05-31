// File: BasisSet/Atom/Gaussian/NR/Imp/Evaluator.C
module;
#include <cassert>
#include <iostream>
#include <blaze/math/DynamicVector.h>

module qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.GaussianIntegrals; 
import qchem.BasisSet.Atom.Evaluators.Gaussian.Internal.ExponentScaler; 
import qchem.Math;

namespace BasisSet::Atom::Evaluators::Gaussian
{
//---------------------------------------------------------------------------
//
//  Start member functions.
//
std::string Evaluator::Name() const
{
    return "SG ";
}

std::string Evaluator::RadialType() const
{
    std::ostringstream os;
    os << "SG";
    return os.str();
}

Cache4*    Evaluator::MakeCache4() const
{
    return new Gaussian_Cache4();
}

rvec_t Evaluator::exponents(size_t N, double emin, double emax, const Irrep_QNs::sym_t& ir)
{
    size_t LMax=3; //TODO how do we get the real LMax(Z) into this?
    ::Gaussian::ExponentScaler ss(N,emin,emax,LMax);
    return ss.Get_es(ir);
}

rvec_t Evaluator::norms() const
{
    size_t N=es.size();    
    rvec_t ret(N);
    for (size_t i=0;i<N;i++) ret[i]=1.0/sqrt(::Gaussian::Integral(2*es[i],2*l)); 
    return ret;
}


rvec_t Evaluator::operator() (const rvec3_t& r) const
{
    return gaussian(norm(r),l,es,ns); 
}

rvec3vec_t Evaluator::Gradient(const rvec3_t& r) const
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

std::ostream&  Evaluator::Write(std::ostream& os) const
{
    return os << " with N=" << es.size() << " basis functions, alpha={" << es[0] << " ... " << es[size()-1] << "}" << std::endl;
}


} //namespace
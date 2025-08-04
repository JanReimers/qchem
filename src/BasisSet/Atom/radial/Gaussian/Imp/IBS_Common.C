// File: Atom/radial/Gaussian/IBS_Common.C  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.
module;
#include <iostream>
#include <iomanip>
#include <iosfwd>
#include <vector>
#include <cassert>
module qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianIntegrals;
import qchem.Symmetry;
import Common.IntPow;
import oml;



namespace Gaussian
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents, Symmetry* sym,size_t l)
    : IrrepBasisSet_Common<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l);
};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents, Symmetry* sym,size_t l, const std::vector<int>& ml)
    : IrrepBasisSet_Common<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,l),l,ml);
};

Vector<double> IrrepBasisSet::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Gaussian::Norm(e,l);
    return ns;
}

IrrepBasisSet::Vec     IrrepBasisSet::operator() (const RVec3& r) const
{
    double mr=norm(r);
    return uintpow(mr,l)*DirectMultiply(ns,exp(-mr*mr*es));
}
IrrepBasisSet::Vec3Vec IrrepBasisSet::Gradient   (const RVec3& r) const
{
    Vec3Vec ret(size());
    double mr=norm(r);
    if (mr==0.0) 
    {
        
        Fill(ret,RVec3(0,0,0));
        return ret; //Cusp at the origin so grad is undefined.
    }
    assert(mr>0);
    Fill(ret,r/mr);
    Vec gr=DirectMultiply(operator()(r),(l/mr-2*mr*es));
    size_t i=0;
    for (auto& ir:ret) ir*=gr(++i);
    return ret;

}


std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << "Spherical Gaussian L=" << *GetSymmetry()
    << " with " << GetNumFunctions() << " basis functions, alpha={";
    // for (auto b:*this) os << *b;
    os << es(1) << " ... " << es(GetNumFunctions());
    os << "}" << std::endl;
    return os;
}




} //namespace

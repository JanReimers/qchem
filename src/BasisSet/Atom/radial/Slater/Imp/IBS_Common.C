// File: Atom/radial/Slater/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Slater functions.
module;
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
module qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import qchem.Symmetry;
import Common.IntPow;
import oml;

namespace Slater
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,Symmetry* sym,size_t L)
    : IrrepBasisSet_Common<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L);

};

IrrepBasisSet::IrrepBasisSet(const Vector<double>& exponents,Symmetry* sym, size_t L, const std::vector<int>& ml)
    : IrrepBasisSet_Common<double>(sym)
    , AtomIrrepIEClient(exponents.size())
{
    Init(exponents,Norms(exponents,L),L,ml);
};

Vector<double> IrrepBasisSet::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Slater::Norm(e,l+1);
    return ns;
}

IrrepBasisSet::Vec     IrrepBasisSet::operator() (const RVec3& r) const
{
    double mr=norm(r);
    return uintpow(mr,l)*DirectMultiply(ns,exp(-mr*es));
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
    Vec gr=DirectMultiply(operator()(r),(l/mr-es));
    size_t i=0;
    for (auto& ir:ret) ir*=gr(++i);
    return ret;

}

std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << "Spherical Slater L=" << *GetSymmetry()
    << " with " << GetNumFunctions() << " basis functions, alpha={";
    // for (auto b:*this) os << *b;
    os << es(1) << " ... " << es(GetNumFunctions());
    os << "}" << std::endl;
    return os;
}

} //namespace 
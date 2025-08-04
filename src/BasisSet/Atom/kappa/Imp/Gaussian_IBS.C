// File: Atom/kappa/Gaussian_IBS.C  Restricted Kinetic Balance (RKB) Irrep Basis Set (IBS).
module;
#include <memory>
#include <iostream>
#include <cassert>
#include <cmath>
module qchem.BasisSet.Atom.Internal.kappa.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianIntegrals;
import qchem.BasisSet.Atom.Internal.radial.GaussianRk;
import qchem.Symmetry.Okmj;
import qchem.Streamable;
import qchem.Orbital_DHF_IBS;
import Common.IntPow;

using std::endl;

namespace Atom_kappa
{
namespace Gaussian
{
  
Orbital_RKB_IBS::Orbital_RKB_IBS
    (const DB_cache<double>* db
        , const Vector<double>& exponents
        , int kappa)
    : Orbital_RKB_IBS_Common<double>
        (db
        , new Omega_k_Sym(kappa)
        , kappa
        , new Large_Orbital_IBS<double>(db,exponents, kappa)
        , new Small_Orbital_IBS<double>(db,exponents,kappa)
        )
{
    
};

std::ostream&  Orbital_RKB_IBS::Write(std::ostream& os) const
{
    os << "Dirac basis set." << endl << "    Large: " << *itsRKBL << endl << "    Small: " << *itsRKBS << endl;
    return os;
}


template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const DB_cache<T>* db,
    const Vector<T>& exponents,int kappa)
    : Orbital_RKBL_IBS_Common<T>(new Omega_k_Sym(kappa),kappa)
    , Orbital_RKBL_IE<T>(db)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    AtomIrrepIEClient::Init(exponents,Norms(exponents,l),l);    
   
};
template <class T> Vector<double> Large_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Gaussian::Norm(e,l);
    return ns;
}
template <class T> Large_Orbital_IBS<T>::Vec     Large_Orbital_IBS<T>::operator() (const RVec3& r) const
{
    double mr=norm(r);
    return uintpow(mr,l)*DirectMultiply(ns,exp(-mr*mr*es));
}

template <class T> Large_Orbital_IBS<T>::Vec3Vec Large_Orbital_IBS<T>::Gradient   (const RVec3& r) const
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

template <class T> std::ostream&  Large_Orbital_IBS<T>::Write(std::ostream& os) const
{
    os << "Gaussian     " << this->GetSymmetry()
    << "               r^" << l << "*exp(-e*r^2), e={";
    for (auto e:es) os << e << ",";
    os << "}";
    return os;
}

template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS(const DB_cache<T>* db,
    const Vector<T>& exponents,int kappa)
    : Orbital_RKBS_IBS_Common<T>(new Omega_k_Sym(-kappa),kappa)
    , Orbital_RKBS_IE<T>(db)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    AtomIrrepIEClient::Init(exponents,Norms(exponents,l),l);  
    
}

template <class T> Vector<double> Small_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=1.0/sqrt(::Gaussian::IE_Primatives::Grad2(e,e,l,l));
    return ns;
}
template <class T> Small_Orbital_IBS<T>::Vec     Small_Orbital_IBS<T>::operator() (const RVec3& r) const
{
    const Large_Orbital_IBS<T>* l1=dynamic_cast<const Large_Orbital_IBS<T>*>(large);
    double mr=norm(r);
    Vec f=-2*es*mr;
    if (l1->kappa >0) 
        f+=(2*l1->kappa+1)/mr;
        
    Vec n=DirectDivide(ns,l1->ns); //Pr(r) is already normalized.
    Vec nf=DirectMultiply(n,f);
    return DirectMultiply(nf,(*large)(r)); 
}

template <class T> Small_Orbital_IBS<T>::Vec3Vec Small_Orbital_IBS<T>::Gradient   (const RVec3& r) const
{
    assert(false);
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
template <class T> std::ostream&  Small_Orbital_IBS<T>::Write(std::ostream& os) const
{
    os << "Gaussian     " << this->GetSymmetry()
    << "               r^" << l << "*exp(-e*r^2), e={";
    for (auto e:es) os << e << ",";
    os << "}";
    return os;
}



}} //namespace

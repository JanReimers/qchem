// File: Atom/kappa/Slater_IBS.C  Slater Irrep Basis Set (IBS) with Restricted Kinetic Balance (RKB).
module;
#include <iostream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <memory>
module qchem.BasisSet.Atom.Internal.kappa.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.ExponentScaler;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import qchem.Symmetry.Okmj;
import qchem.Streamable;
import qchem.DHF_IBS;
import Common.IntPow;

using std::endl;

namespace Atom_kappa
{
namespace Slater
{
//
//  Concrete  Slater basis set.
//

Orbital_IBS::Orbital_IBS
    (const DB_cache<double>* db
        , const Vector<double>& exponents
        , int kappa)
    : Orbital_RKB_IBS_Common1<double>
        (db
        , new Omega_k_Sym(kappa)
        , kappa
        , new Large_Orbital_IBS<double>(db,exponents, kappa)
        , new Small_Orbital_IBS<double>(db,exponents,kappa)
        )
{
    
};


std::ostream&  Orbital_IBS::Write(std::ostream& os) const
{
    os << "Dirac basis set." << endl << "    Large: " << *itsRKBL << endl << "    Small: " << *itsRKBS << endl;
    return os;
}


//-----------------------------------------------------------------------------------------------
//
//  Large sector
//
template <class T> Large_Orbital_IBS<T>::Large_Orbital_IBS(const DB_cache<T>* db,
        const Vector<T>& exponents,int kappa)
    : Orbital_RKBL_IBS_Common1<T>(new Omega_k_Sym(kappa),kappa)
    , Orbital_RKBL_IE<T>(db)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    Init(exponents,Norms(exponents,l),l);
   

};

template <class T> Vector<double> Large_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=::Slater::Norm(e,l+1);
    return ns;
}

template <class T> Large_Orbital_IBS<T>::Vec     Large_Orbital_IBS<T>::operator() (const RVec3& r) const
{
    double mr=norm(r);
    return uintpow(mr,l)*DirectMultiply(ns,exp(-mr*es));
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
    Vec gr=DirectMultiply(operator()(r),(l/mr-es));
    size_t i=0;
    for (auto& ir:ret) ir*=gr(++i);
    return ret;
}

template <class T> std::ostream&  Large_Orbital_IBS<T>::Write(std::ostream& os) const
{
    os << "Slater     " << this->GetSymmetry()
    << "             r^" << l << "*exp(-e*r), e={";
    for (auto e:es) os << e << ",";
    os << "}";
    return os;
}

//-----------------------------------------------------------------------------------------------
//
//  Small sector
//
template <class T> Small_Orbital_IBS<T>::Small_Orbital_IBS
    (const DB_cache<double>* db,
        const Vector<T>& exponents
        , int kappa)
    : Orbital_RKBS_IBS_Common1<T>(new Omega_k_Sym(-kappa),kappa)
    , Orbital_RKBS_IE<T>(db)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    Init(exponents,Norms(exponents,l),l);
};

template <class T> Vector<double> Small_Orbital_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    for (auto e:es) ns(++i)=1.0/sqrt(::Slater::IE_Primatives::Grad2(e,e,l,l));
    return ns;
}
template <class T> Small_Orbital_IBS<T>::Vec     Small_Orbital_IBS<T>::operator() (const RVec3& r) const
{
    const Large_Orbital_IBS<T>* l1=dynamic_cast<const Large_Orbital_IBS<T>*>(large);
    double mr=norm(r);
    Vec f=-es;
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
    Vec gr=DirectMultiply(operator()(r),(l/mr-es));
    size_t i=0;
    for (auto& ir:ret) ir*=gr(++i);
    return ret;
}
template <class T> std::ostream&  Small_Orbital_IBS<T>::Write(std::ostream& os) const
{
    os << "Slater RKB " << this->GetSymmetry();
    if (kappa>0)
        os << "[ " << std::setw(2) << 2*kappa+1 << "/r - e ]";
    else
        os << "[       -e ]";
    os << "*r^" << l << "*exp(-e*r), e={";
    for (auto e:es) os << e << ",";
    os << "}";
    return os;
}


}} //namespace

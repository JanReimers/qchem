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
import BasisSet.Atom.Gaussian_IBS;
import qchem.Symmetry.Okmj;
import qchem.Streamable;
import qchem.Orbital_DHF_IBS;
import Common.IntPow;

using std::endl;

namespace Atom_kappa
{
namespace Gaussian
{

double IE_Primatives_sgkappa::Inv_r1(double ea, double eb,size_t l_total) const
{
    return 4*ea*eb*::Gaussian::Integral(ea+eb,l_total+1); //Don't count the r^2 in dr^3

}

Orbital_RKB_IBS::Orbital_RKB_IBS
    (const DB_cache<double>* db,const ::IE_Primatives* pie, const Vector<double>& exponents, int kappa)
    : IrrepBasisSet_Common<double>(new Omega_k_Sym(kappa))
    , Orbital_RKB_IBS_Common<double>
        (db, kappa
        , new Orbital_RKBL_IBS<double>(db,pie ,new Gaussian_IBS(exponents,Omega_k_Sym::l(kappa),{}),exponents, kappa)
        , new Orbital_RKBS_IBS<double>(db,new Gaussian_RKBS_IBS(exponents,Omega_k_Sym::l(kappa),{}),exponents,kappa) //Known memory leak, need redesign
        )
{
    
};

std::ostream&  Orbital_RKB_IBS::Write(std::ostream& os) const
{
    os << "Dirac basis set." << endl << "    Large: " << *itsRKBL << endl << "    Small: " << *itsRKBS << endl;
    return os;
}

//-----------------------------------------------------------------------------------------------
//
//  Large sector
//
template <class T> Orbital_RKBL_IBS<T>::Orbital_RKBL_IBS
(const DB_cache<T>* db,const IE_Primatives* pie,const IBS_Evaluator* eval,const Vector<T>& exponents,int kappa)
    : ::Gaussian::IrrepBasisSet(exponents,new Omega_k_Sym(kappa),Omega_kmj_Sym::l(kappa))
    , Atom::Orbital_RKBL_IBS<T>(db,pie,eval,kappa)
{};



//-----------------------------------------------------------------------------------------------
//
//  Small sector
//
template <class T> Orbital_RKBS_IBS<T>::Orbital_RKBS_IBS(const DB_cache<T>* db,const IBS_Evaluator* eval,
    const Vector<T>& exponents,int kappa)
    : IrrepBasisSet_Common<T> (new Omega_k_Sym(-kappa))
    , Atom::Orbital_RKBS_IBS<T>(db,eval,kappa)
    , AtomIrrepIEClient(exponents.size())
{
    size_t l=Omega_kmj_Sym::l(kappa);
    AtomIrrepIEClient::Init(exponents,Norms(exponents,l),l);  
    
}

template <class T> Vector<double> Orbital_RKBS_IBS<T>::Norms(const Vector<double>& es, size_t l) const
{
    Vector<double> ns(es.size());
    int i=0;
    ::Gaussian::IE_Primatives pie;
    for (auto e:es) ns(++i)=1.0/sqrt(pie.Grad2(e,e,l,l));
    return ns;
}
template <class T> Orbital_RKBS_IBS<T>::Vec     Orbital_RKBS_IBS<T>::operator() (const RVec3& r) const
{
    const Orbital_RKBL_IBS<T>* l1=dynamic_cast<const Orbital_RKBL_IBS<T>*>(large);
    double mr=norm(r);
    Vec f=-2*es*mr;
    if (l1->kappa >0) 
        f+=(2*l1->kappa+1)/mr;
        
    Vec n=DirectDivide(ns,l1->ns); //Pr(r) is already normalized.
    Vec nf=DirectMultiply(n,f);
    return DirectMultiply(nf,(*large)(r)); 
}

template <class T> Orbital_RKBS_IBS<T>::Vec3Vec Orbital_RKBS_IBS<T>::Gradient   (const RVec3& r) const
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
template <class T> std::ostream&  Orbital_RKBS_IBS<T>::Write(std::ostream& os) const
{
    os << "Gaussian     " << this->GetSymmetry()
    << "               r^" << l << "*exp(-e*r^2), e={";
    for (auto e:es) os << e << ",";
    os << "}";
    return os;
}



}} //namespace

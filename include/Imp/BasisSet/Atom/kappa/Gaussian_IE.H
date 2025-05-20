// File: Atom/kappa/Gaussian_IE.H  Integral Engine for RKB Gaussians.
#ifndef _Atom_kappa_Gaussian_IE_H_
#define _Atom_kappa_Gaussian_IE_H_


#include "../IE.H"
#include "Imp/BasisSet/Atom/l/Gaussian_IE.H"

class DiracIntegralTests;

namespace Atom_kappa
{
namespace Gaussian
{
    
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public virtual ::Gaussian::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
    , public virtual ::Gaussian::IE_Primatives
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db) : AtomIE_RKBS<T>(db) {};
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
};


}} //namespace SphericalGaussian_RKB


#endif //_Atom_kappa_Gaussian_IE_H_

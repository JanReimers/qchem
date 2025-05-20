// File: Atom/kappa/Slater_IE.H  Inegral Engine (IE) for Slater basis set with Restricted Kinetic Balance (RKB).
#ifndef _Atom_kappa_Slater_IE_H_
#define _Atom_kappa_Slater_IE_H_

#include "../radial//Slater/IE_Primatives.H"
#include "../IE_DHF.H"

class DiracIntegralTests;

namespace Atom_kappa
{
namespace Slater
{
    
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public virtual ::Slater::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
    , public virtual ::Slater::IE_Primatives
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db) : AtomIE_RKBS<T>(db) {};
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
};




}} //namespace

#endif //_Atom_kappa_Slater_IE_H_

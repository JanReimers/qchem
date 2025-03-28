// File: Slater_mj/IntegralEngine.H  Integral Engine for Slater function with kappa,mj QNs.
#ifndef _Slater_mj_IntegralEngine_H_
#define _Slater_mj_IntegralEngine_H_

#include "../Atom/radial//Slater/IE_Primatives.H"
#include "../Atom_IE_DHF.H"

class DiracIntegralTests;

namespace Slater_mj
{
    
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public virtual Slater::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
    , public virtual Slater::IE_Primatives
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db) : AtomIE_RKBS<T>(db) {};
    virtual double Grad2(double ea , double eb,size_t la,size_t lb) const;
    virtual double Nuclear(double ea , double eb,size_t l_total) const;
};




} //namespace

#endif //_Slater_m_IntegralEngine_H_

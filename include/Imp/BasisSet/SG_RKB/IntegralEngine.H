// File: SphericalGaussian_RKB_IE.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussian_RKB_IE1_H_
#define _SphericalGaussian_RKB_IE1_H_


#include "IEClient.H"
#include "../AtomIE.H"
#include "../SphericalGaussian/IntegralEngine.H"

class DiracIntegralTests;

namespace SphericalGaussian_RKB
{
    
template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public SphericalGaussian::IE_Primatives
{
protected:
    Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
};

template <class T> class Orbital_RKBS_IE
    : public     AtomIE_RKBS<T>
{
protected:
    Orbital_RKBS_IE(const DB_cache<double>* db) : AtomIE_RKBS<T>(db) {};
    virtual double Kinetic(double ea , double eb,size_t la,size_t lb) const;
    virtual double Nuclear(double ea , double eb,size_t l_total) const;
};


} //namespace SphericalGaussian_RKB


#endif //_SphericalGaussian_RKB_IE_H_

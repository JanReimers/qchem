// File: SphericalGaussian_RKB_IE.H  Integral Engine for Spherical gaussians.
#ifndef _SphericalGaussian_RKB_IE1_H_
#define _SphericalGaussian_RKB_IE1_H_


#include <AnalyticIE.H>
#include "IEClient.H"
#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/SphericalGaussian/IntegralEngine.H"
#include "oml/vector.h"

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


} //namespace SphericalGaussian_RKB


#endif //_SphericalGaussian_RKB_IE_H_

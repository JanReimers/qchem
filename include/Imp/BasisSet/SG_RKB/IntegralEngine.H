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
    : public virtual AtomIE_1EFit_Common
    , public     AtomIE_RKBL<T>
    , public SphericalGaussian::IE_Primatives
    {
    protected:
        Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
    };

    
// class Small_IntegralEngine1
//     : public virtual Integrals_1E<double>
//     , public AtomIE_1E<double>
//     , public DB_RKBS<double>
// {
// public:
//     virtual  double Integral(qchem::IType,double ea , double eb,size_t l) const;
// };



} //namespace SphericalGaussian_RKB


#endif //_SphericalGaussian_RKB_IE_H_

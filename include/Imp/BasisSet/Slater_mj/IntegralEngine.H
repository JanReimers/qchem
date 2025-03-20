// File: Slater_mj/IntegralEngine.H  Integral Engine for Slater function with kappa,mj QNs.
#ifndef _Slater_mj_IntegralEngine_H_
#define _Slater_mj_IntegralEngine_H_

#include "Imp/BasisSet/Dirac_IBS.H"
#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/Slater/IntegralEngine.H"

class DiracIntegralTests;

namespace Slater_mj
{
    
    template <class T> class Orbital_RKBL_IE
    : public     AtomIE_RKBL<T>
    , public Slater::IE_Primatives
    {
    protected:
        Orbital_RKBL_IE(const DB_cache<double>* db) : AtomIE_RKBL<T>(db) {};
    };



} //namespace

#endif //_Slater_m_IntegralEngine_H_

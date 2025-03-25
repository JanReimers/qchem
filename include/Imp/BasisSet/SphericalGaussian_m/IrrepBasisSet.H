// File: SphericalGaussian_m/IrrepBasisSet.H  Spherical Slater basis set with orbital angular momentum l,m.
#ifndef _SphericalGaussian_m_IrrepBasisSet_H_
#define _SphericalGaussian_m_IrrepBasisSet_H_



#include "IEClient.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../SphericalGaussian/IntegralEngine.H"
#include "../HF_IBS_Common.H"
#include "oml/vector.h"
#include <set>

namespace SphericalGaussian_m
{

    // Common implementation for orbital and fit basis sets.
class IrrepBasisSet
: public virtual ::IrrepBasisSet
, public         IrrepBasisSetCommon
, public         IrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t L, int m);
    IrrepBasisSet(const std::set<double>& exponents, size_t L,int m);
    virtual std::ostream&  Write(std::ostream&    ) const;

};

class Orbital_IBS
: public virtual TOrbital_HF_IBS<double>
// , public virtual TOrbital_DFT_IBS<double>
, public         IrrepBasisSet
, public         Orbital_IBS_Common<double>
// , public         Orbital_DFT_IBS_Common<double>
, public         Orbital_HF_IBS_Common<double>
, public         SphericalGaussian::Orbital_IE

{
public:
Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, int m)
    : IrrepBasisSet(exponents,L,m)
    , Orbital_IBS_Common<double>()
    , SphericalGaussian::Orbital_IE(db)
{};

Orbital_IBS(const DB_BS_2E<double>* db,const std::set<double>& exponents, size_t L, int m)
    : IrrepBasisSet(exponents,L,m)
    , Orbital_IBS_Common<double>()
    , SphericalGaussian::Orbital_IE(db)
{};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};


} //namespace SphericalGaussian_m

#endif //_SphericalGaussian_m_IrrepBasisSet_H_

// File: Slater_m/IrrepBasisSet.H  Spherical Slater basis set with orbital angular momentum l,m.
#ifndef _Slater_m_IrrepBasisSet_H_
#define _Slater_m_IrrepBasisSet_H_



#include "IEClient.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../DFT_IBS_Common.H"
#include "../HF_IBS_Common.H"
#include "../Fit_IBS_Common.H"
#include "../Atom/l/Slater_IE.H"
#include "oml/vector.h"

namespace Slater_m
{
// Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         IrrepBasisSetCommon
    , public         IrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t L, int m);
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
};

class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public virtual TOrbital_DFT_IBS<double>
    , public         IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Atoml::Slater::Orbital_IE

{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, int m)
    : IrrepBasisSet(exponents,L,m)
    , Orbital_IBS_Common<double>()
    , Atoml::Slater::Orbital_IE(db)
    {};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};
 

} //namespace Slater_m

#endif //_Slater_m_IrrepBasisSet_H_

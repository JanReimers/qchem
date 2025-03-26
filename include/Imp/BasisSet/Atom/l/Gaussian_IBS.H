// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.
#ifndef _Atom_l_Gaussian_IBS_H_
#define _Atom_l_Gaussian_IBS_H_



#include "Gaussian_IEC.H"
#include "GaussianIE.H"
#include "../../IrrepCommon.H"
#include "../../TIrrepCommon.H"
#include "../../DFT_IBS_Common.H"
#include "../../HF_IBS_Common.H"
#include "../../Fit_IBS_Common.H"
#include "oml/vector.h"

namespace Atoml
{
namespace Gaussian
{
   
// Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         IrrepBasisSetCommon
    , public         IrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t L);
    virtual std::ostream&  Write(std::ostream&    ) const;
};

class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Orbital_IE
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , Orbital_IBS_Common<double>()
    , Orbital_IE(db)
    {};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};

class Fit_IBS 
: public virtual ::Fit_IBS
, public         IrrepBasisSet
, public         TIrrepBasisSetCommon<double>
, public Fit_IBS_Common
, public Fit_IE
{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , TIrrepBasisSetCommon<double>()
    , Fit_IE(db)
    {};
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
};


} //namespace
} //namespace

#endif //_Atom_l_Gaussian_IBS_H_

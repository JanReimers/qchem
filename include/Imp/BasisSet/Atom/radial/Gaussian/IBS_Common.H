// File: Atom/radial/Gaussian/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.
#ifndef _Atom_radial_Gaussian_IBS_Common_H_
#define _Atom_radial_Gaussian_IBS_Common_H_
#include "../../../AtomIEClient.H"
#include "../../../IrrepCommon.H"
#include "../../../TIrrepCommon.H"
#include "../../../DFT_IBS_Common.H"
#include "../../../HF_IBS_Common.H"
#include "../../../Fit_IBS_Common.H"
#include "../../l/GaussianIE.H"

#include "oml/vector.h"

namespace Gaussian
{
    // Common implementation for orbital and fit basis sets.
class IrrepBasisSet
: public virtual ::IrrepBasisSet
, public         IrrepBasisSetCommon
, public         AtomIrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t l);
    IrrepBasisSet(const Vector<double>& exponents, size_t l, int m);
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
};

class Fit_IBS 
: public virtual ::Fit_IBS
, public         IrrepBasisSet
, public         TIrrepBasisSetCommon<double>
, public Fit_IBS_Common
, public Atoml::Gaussian::Fit_IE
{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , TIrrepBasisSetCommon<double>()
    , Atoml::Gaussian::Fit_IE(db)
    {};
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
};


} //namespace

#endif //_Atom_radial_Gaussian_IBS_Common_H_

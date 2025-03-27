// File: Atom/ml/Gaussian_IBS.H  r^l exp(-ar^2)*Y_lm type Irrep Basis set (IBS).
#ifndef _Atom_ml_Gaussian_IBS_H_
#define _Atom_ml_Gaussian_IBS_H_


#include "../../AtomIEClient.H"
#include "../../IrrepCommon.H"
#include "../../TIrrepCommon.H"
#include "../l/GaussianIE.H"
#include "../../HF_IBS_Common.H"
#include "oml/vector.h"
#include <set>

namespace Atom_ml
{
namespace Gaussian
{

    // Common implementation for orbital and fit basis sets.
class IrrepBasisSet
: public virtual ::IrrepBasisSet
, public         IrrepBasisSetCommon
, public         AtomIrrepIEClient
{
public:
    IrrepBasisSet(const Vector<double>& exponents, size_t L, int m);
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
};

class Orbital_IBS
: public virtual TOrbital_HF_IBS<double>
// , public virtual TOrbital_DFT_IBS<double>
, public         IrrepBasisSet
, public         Orbital_IBS_Common<double>
// , public         Orbital_DFT_IBS_Common<double>
, public         Orbital_HF_IBS_Common<double>
, public         Atoml::Gaussian::Orbital_IE

{
public:
Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, int m)
    : IrrepBasisSet(exponents,L,m)
    , Orbital_IBS_Common<double>()
    , Atoml::Gaussian::Orbital_IE(db)
{};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};


} //namespace 
} //namespace 

#endif //_Atom_ml_Gaussian_IBS_H_

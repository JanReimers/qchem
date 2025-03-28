// File: Atom/ml/Gaussian_IBS.H  r^l exp(-ar^2)*Y_lm type Irrep Basis set (IBS).
#ifndef _Atom_ml_Gaussian_IBS_H_
#define _Atom_ml_Gaussian_IBS_H_


#include "../radial/Gaussian/IBS_Common.H"
#include "../l/GaussianIE.H"
#include "../../HF_IBS_Common.H"
#include "../../TIrrepCommon.H"

namespace Atom_ml
{
namespace Gaussian
{

class Orbital_IBS
: public virtual TOrbital_HF_IBS<double>
// , public virtual TOrbital_DFT_IBS<double>
, public         ::Gaussian::IrrepBasisSet
, public         Orbital_IBS_Common<double>
// , public         Orbital_DFT_IBS_Common<double>
, public         Orbital_HF_IBS_Common<double>
, public         Atoml::Gaussian::Orbital_IE

{
public:
Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, int m)
    : ::Gaussian::IrrepBasisSet(exponents,L,m)
    , Orbital_IBS_Common<double>()
    , Atoml::Gaussian::Orbital_IE(db)
{
    InsertBasisFunctions();   

};

    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};


} //namespace 
} //namespace 

#endif //_Atom_ml_Gaussian_IBS_H_

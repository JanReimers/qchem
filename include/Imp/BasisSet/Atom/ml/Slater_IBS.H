// File: Atom/ml/Slater_IBS.H  Slater Irrep Basis Set (IBS) with orbital angular momentum l,m.
#ifndef _Atom_ml_Slater_IBS_H_
#define _Atom_ml_Slater_IBS_H_



#include "../../AtomIEClient.H"
#include "../../IrrepCommon.H"
#include "../../TIrrepCommon.H"
#include "../../DFT_IBS_Common.H"
#include "../../HF_IBS_Common.H"
#include "../../Fit_IBS_Common.H"
#include "../l/Slater_IE.H"
#include "oml/vector.h"

namespace Atom_ml
{
namespace Slater
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
 

}} //namespace Atom_ml::Slater

#endif //_Atom_ml_Slater_IBS_H_

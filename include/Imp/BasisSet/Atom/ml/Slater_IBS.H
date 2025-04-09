// File: Atom/ml/Slater_IBS.H  Slater Irrep Basis Set (IBS) with orbital angular momentum l,m.
#ifndef _Atom_ml_Slater_IBS_H_
#define _Atom_ml_Slater_IBS_H_

#include "../l/Slater_IE.H"
#include "../radial/Slater/IBS_Common.H"
#include "../../TIBS_Common.H"
#include "../../DFT_IBS_Common.H"
#include "../../HF_IBS_Common.H"
#include "../../Fit_IBS_Common.H"

namespace Atom_ml
{
namespace Slater
{

class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public virtual TOrbital_DFT_IBS<double>
    , public         ::Slater::IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Atoml::Slater::Orbital_IE

{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, int m);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();

};
 

}} //namespace Atom_ml::Slater

#endif //_Atom_ml_Slater_IBS_H_

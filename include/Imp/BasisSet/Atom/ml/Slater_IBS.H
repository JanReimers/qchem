// File: Atom/ml/Slater_IBS.H  Slater Irrep Basis Set (IBS) with orbital angular momentum l,m.
#ifndef _Atom_ml_Slater_IBS_H_
#define _Atom_ml_Slater_IBS_H_

#include "../radial/Slater/IBS_Common.H"
import qchem.BasisSet.Atom.l.SlaterBS;
import qchem.BasisSet.qchem.BasisSet.IBS_Common;
import qchem.BasisSet.Imp.HeapDB;
import qchem.HF_IBS;
import qchem.DFT_IBS;

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
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const ::IrrepBasisSet::RVec3&) const;
private:
    void InsertBasisFunctions();

};
 

}} //namespace Atom_ml::Slater

#endif //_Atom_ml_Slater_IBS_H_

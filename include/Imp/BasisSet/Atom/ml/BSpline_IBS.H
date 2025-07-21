// File: Atom/ml/BSpline_IBS.H  BSpline Irrep Basis Set (IBS) with orbital angular momentum l,m.
#ifndef _Atom_ml_BSpline_IBS_H_
#define _Atom_ml_BSpline_IBS_H_

#include "../radial/BSpline/IBS_Common.H"

import qchem.BasisSet.Atom.l.BSplineBS;
import qchem.BasisSet.qchem.BasisSet.IBS_Common;
import qchem.BasisSet.Atom.l.GaussianBS;

namespace Atom_ml
{
namespace BSpline
{

template <size_t K> class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    // , public virtual TOrbital_DFT_IBS<double>
    , public         ::BSpline::IrrepBasisSet<K>
    , public         Orbital_IBS_Common<double>
    // , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Atoml::BSpline::Orbital_IE<K>

{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L,  const std::vector<int>& ml);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const ::IrrepBasisSet::RVec3&) const;
private:
    void InsertBasisFunctions();

};
 

}} //namespace Atom_ml::BSpline

#endif //_Atom_ml_BSpline_IBS_H_

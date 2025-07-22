// File: Atom/ml/BSpline_BS.H BSpline Basis Set for atoms, no m degeneracy.
module;
#include <bspline/Core.h>
#include <iosfwd>

export module qchem.BasisSet.Atom.ml.BSplineBS;
import qchem.BasisSet.Atom.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.BS_Common;

import qchem.BasisSet.IBS_Common;
import qchem.BasisSet.Common;

import qchem.BasisFunction;
import qchem.BasisSet.Atom.ml.Angular;
import qchem.HF_IBS;
import qchem.BasisSet.Imp.HeapDB;
import qchem.Fit_IBS;
import qchem.BasisSet;


export namespace Atom_ml
{
namespace BSpline
{
template <size_t K> class BasisFunction
    : public Atoml::BSpline::BasisFunction<K>
{
    typedef Atoml::BSpline::BasisFunction<K> Base;
    typedef typename Base::spline_t spline_t;
public:
    // BasisFunction();
    BasisFunction(const spline_t&, int l, int ml, double norm);
    
    virtual ::Real_BF* Clone(        ) const;

private:
    int      ml;
};
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
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();

};
template <size_t K> class BasisSet 
    : public ::BSpline::BS_Common<K>
    , public IE_BS_2E_Angular
{
public:
    BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec);
    
};

}} //namespace Atom_ml::BSpline

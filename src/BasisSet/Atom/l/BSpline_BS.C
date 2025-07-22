// File: Atom/l/BSpline_BS.H BSpline Basis Set for atoms.
module;
#include <bspline/Core.h>
#include <iosfwd>

export module qchem.BasisSet.Atom.l.BSplineBS;
import qchem.BasisSet;
import qchem.HF_IBS;
import qchem.Fit_IBS;
import qchem.BasisSet.IBS_Common;
import qchem.BasisSet.Common;
import qchem.BasisFunction;
import qchem.BasisSet.Atom.l.Angular;
import qchem.BasisSet.Imp.HeapDB;

import qchem.Basisset.Atom.radial.BSpline.BS_Common;
import qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.Basisset.Atom.radial.BSpline.IE_Primatives;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;

export namespace Atoml
{
namespace BSpline
{
    // Basis function
template <size_t K> class BasisFunction
    : public virtual ::Real_BF
{
protected:
    typedef bspline::Spline<double, K> spline_t;
public:
    // BasisFunction();
    BasisFunction(const spline_t&, int l, double norm);
    
    virtual std::ostream&    Write(std::ostream&) const;
    virtual ::Real_BF* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    spline_t itsSpline;
    bspline::Spline<double, K-1> itsDxSpline;
    int      itsL;
    double   itsNormalization;
};

    // Integral engine
template <size_t K> class IE_Common
    : public virtual ::BSpline::IE_Primatives<K>
    , public ::BSpline::IE_Overlap<double,K>
    , public ::BSpline::IE_Kinetic  <double,K>
    , public ::BSpline::IE_Inv_r1<double,K>
{
public:
protected:
    IE_Common(const DB_cache<double>* db) : ::BSpline::IE_Overlap<double,K>(db), ::BSpline::IE_Kinetic<double,K>(db), ::BSpline::IE_Inv_r1<double,K>(db) {};
};

template <size_t K> class Orbital_IE
: public IE_Common<K>
, public DB_2E<double>
, public ::BSpline::IE_DFT<double,K>
{
public:
protected:
    Orbital_IE(const DB_BS_2E<double>* db) : IE_Common<K>(db), DB_2E<double>(db), ::BSpline::IE_DFT<double,K>(db) {};

};

template <size_t K> class Fit_IE
: public ::BSpline::IE_Fit<K>
, public ::BSpline::IE_Overlap<double,K>
, public virtual ::BSpline::IE_Primatives<K>

{
protected:
    Fit_IE(const DB_cache<double>* db) : ::BSpline::IE_Fit<K>(db), ::BSpline::IE_Overlap<double,K>(db) {};
};

    // Irrep basis set
template <size_t K> class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         ::BSpline::IrrepBasisSet<K>
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Orbital_IE<K>
{
public:
    using RVec3=::IrrepBasisSet::RVec3;
    Orbital_IBS(const DB_BS_2E<double>* db,size_t N, double rmin, double rmax, size_t L);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

template <size_t K> class Fit_IBS 
: public virtual ::Fit_IBS
, public virtual FitIntegrals
, public ::BSpline::IrrepBasisSet<K>
, public TIBS_Common<double>
, public Fit_IBS_Common
, public Fit_IE<K>

{
public:
    using RVec3=::IrrepBasisSet::RVec3;
    Fit_IBS(const DB_cache<double>* db,size_t N, double rmin, double rmax, size_t L);
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

    // Full basis set.

template <size_t K> class BasisSet 
    : public ::BSpline::BS_Common<K>
    , public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet(size_t N, double rmin, double rmax, size_t Lmax); 
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
   
};

}} //namespace Atoml::BSpline


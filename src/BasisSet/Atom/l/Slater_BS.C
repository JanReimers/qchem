// File: Atom/l/Slater_BS.H Slater Basis Set for atoms.
module;
#include <iosfwd>
#include <vector>

export module qchem.BasisSet.Atom.l.SlaterBS;
import qchem.BasisSet.Atom.radial.SlaterBS;
import qchem.BasisSet.Atom.radial.Slater.IE_Primatives;

import qchem.BasisSet;
import qchem.BasisSet.IBS_Common;
import qchem.HF_IBS;
import qchem.BasisFunction;
import qchem.BasisSet.Atom.l.Angular;
import qchem.BasisSet.Common;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace Slater
{
    // Basis function
class BasisFunction
    : public virtual ::Real_BF
{
public:
    BasisFunction();
    BasisFunction(double ex, int n, int l, double norm);
    
    virtual std::ostream&  Write(std::ostream&) const;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;

private:
    double itsExponent;
    int    itsN;
    int    itsL;
    double itsNormalization;
};

    // Integral engine
class IE_Common
    : public virtual ::Slater::IE_Primatives
    , public AtomIE_Overlap<double>
    , public AtomIE_Kinetic<double>
    , public AtomIE_Nuclear<double>
{
public:
protected:
    IE_Common(const DB_cache<double>* db) : AtomIE_Overlap<double>(db), AtomIE_Kinetic<double>(db), AtomIE_Nuclear<double>(db) {};
};

class Orbital_IE
: public IE_Common
, public DB_2E<double>
, public AtomIE_DFT<double>
{
public:
protected:
    Orbital_IE(const DB_BS_2E<double>* db) : IE_Common(db), DB_2E<double>(db), AtomIE_DFT<double>(db) {};

};

class Fit_IE
: public AtomIE_Fit
, public AtomIE_Overlap<double>
, public virtual ::Slater::IE_Primatives

{
protected:
    Fit_IE(const DB_cache<double>* db) : AtomIE_Fit(db), AtomIE_Overlap<double>(db) {};
};

    // Irrep basis set
class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         ::Slater::IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Orbital_IE
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

class Fit_IBS 
: public virtual ::Fit_IBS
, public virtual FitIntegrals
, public ::Slater::IrrepBasisSet
, public TIBS_Common<double>
, public Fit_IBS_Common
, public Fit_IE

{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L);
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};


    // Full basis set.

class BasisSet 
    : public ::Slater::BS_Common
    , public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax); 
};

}} //namespace Slater::l


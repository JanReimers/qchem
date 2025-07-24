// File: Atom/l/Gaussian_BS.H Gaussian Basis Set for atoms.
module;
#include <iosfwd>
#include <vector>

export module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.IE_Primatives;
import qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.Internal.l.Angular;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.HF_IBS;
import qchem.BasisFunction;

export namespace Atoml
{
namespace Gaussian
{

    // Basis function
class BasisFunction
    : public virtual ::Real_BF
{
public:
    BasisFunction(                             );
    BasisFunction(double theExponent, int theL, double norm);
    virtual ~BasisFunction();
    
    virtual std::ostream&  Write(std::ostream&) const;
    virtual BasisFunction* Clone(        ) const;

    virtual double operator()(const RVec3&) const;
    virtual RVec3  Gradient  (const RVec3&) const;

private:
    double itsExponent;
    int    itsL;
    double itsNormalization;
};

    // Integral engine
class Orbital_IE
: public ::Gaussian::IE_Primatives
, public AtomIE_Overlap<double>
, public AtomIE_Kinetic<double>
, public AtomIE_Nuclear<double>
, public DB_2E<double>
, public AtomIE_DFT<double>
{
protected:
    Orbital_IE(const DB_BS_2E<double>* db) 
        : AtomIE_Overlap<double>(db)
        , AtomIE_Kinetic<double>(db)
        , AtomIE_Nuclear<double>(db)
        , DB_2E<double>(db)
        , AtomIE_DFT<double>(db) {};

};

class Fit_IE
: public AtomIE_Fit
, public AtomIE_Overlap<double>
, public ::Gaussian::IE_Primatives
{
protected:
    Fit_IE(const DB_cache<double>* db) : AtomIE_Fit(db), AtomIE_Overlap<double>(db) {};
};

    // Irrep basis set
class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public         ::Gaussian::IrrepBasisSet
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
, public         ::Gaussian::IrrepBasisSet
, public         TIBS_Common<double>
, public Fit_IBS_Common
, public Atoml::Gaussian::Fit_IE
{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L);
   
    virtual ::Fit_IBS* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};

    // Full basis set.

class BasisSet 
: public ::Gaussian::BS_Common
, public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
};

} //namespace
} //namespace

// File: Atom/l/Slater_BS.H Slater Basis Set for atoms.
module;
#include <iosfwd>
#include <vector>

export module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.IE_Primatives;

import qchem.BasisSet;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.Orbital_HF_IBS;

import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace Slater
{

    // Integral engine
class IE_Common
    : public virtual ::Slater::IE_Primatives
    , public AtomIE_Overlap<double>
    , public AtomIE_Kinetic<double>
    , public AtomIE_Nuclear<double>
{
public:
protected:
    IE_Common(const DB_cache<double>* db,const IE_Primatives* pie) 
    : AtomIE_Overlap<double>(db,pie)
    , AtomIE_Kinetic<double>(db,pie)
    , AtomIE_Nuclear<double>(db,pie) {};
};


class Fit_IE
: public AtomIE_Fit
, public AtomIE_Overlap<double>
, public ::Slater::IE_Primatives

{
protected:
    Fit_IE(const DB_cache<double>* db,const IE_Primatives* pie) 
    : AtomIE_Fit(db)
    , AtomIE_Overlap<double>(db,pie) {};
};

    // Irrep basis set
class Orbital_IBS
    : public virtual Orbital_HF_IBS<double>
    , public         ::Slater::IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         IE_Common
    , public AtomIE_DFT<double>
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L, const std::vector<int>& ml);
    
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;

};

class Fit_IBS 
: public virtual ::Fit_IBS
, public virtual FitIntegrals
, public ::Slater::IrrepBasisSet
, public Fit_IBS_Common
, public Fit_IE

{
public:
    Fit_IBS(const DB_cache<double>* db,const IE_Primatives* pie,const Vector<double>& exponents, size_t L);
   
};


    // Full basis set.

class BasisSet 
    : public ::Slater::BS_Common
    , public ::Slater::IE_Primatives
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax); 
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
};

}} //namespace Slater::l


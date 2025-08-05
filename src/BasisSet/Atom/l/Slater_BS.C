// File: Atom/l/Slater_BS.H Slater Basis Set for atoms.
module;
#include <iosfwd>
#include <vector>

export module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.IE_Primatives;
import qchem.BasisSet.Atom.Orbital_1E_IBS;
import qchem.BasisSet.Atom.Orbital_DFT_IBS;
import qchem.BasisSet.Atom.Fit_IBS;
import qchem.BasisSet;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.Orbital_HF_IBS;

import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace Slater
{

// Irrep basis set
class Orbital_IBS
    : public virtual Orbital_HF_IBS<double>
    , public         ::Slater::IrrepBasisSet
    , public         Orbital_IBS_Common<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public Atom::Orbital_IBS<double>
    , public Atom::Orbital_DFT_IBS<double>
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
, public Atom::Fit_IBS

{
public:
    Fit_IBS(const DB_cache<double>* db,const ::IE_Primatives* pie,const Vector<double>& exponents, size_t L);
   
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


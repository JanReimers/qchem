// File: BasisSet/Atom/l/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.IE_Primatives;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.IBS;

export namespace Atoml
{
namespace Slater
{
class Orbital_IBS
    : public ::Slater::IrrepBasisSet
    , public Atom::Orbital_HF_IBS <double>
    , public Atom::Orbital_IBS    <double>
    , public Atom::Orbital_DFT_IBS<double>
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const ::IE_Primatives* pie,const Vector<double>& exponents, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,const ::IE_Primatives* pie,const Vector<double>& exponents, size_t L, const std::vector<int>& ml);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
};

class Fit_IBS 
: public ::Slater::IrrepBasisSet
, public Atom::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,const ::IE_Primatives* pie,const Vector<double>& exponents, size_t L);
};

class BasisSet 
    : public ::Slater::BS_Common
    , public ::Slater::IE_Primatives
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax); 
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
};

}} //namespace Slater::l


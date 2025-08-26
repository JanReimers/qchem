// File: Atom/l/Gaussian_BS.H Gaussian Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import BasisSet.Atom.Gaussian_IBS;
import BasisSet.Atom.Gaussian_BS;
import BasisSet.Atom.IBS_Evaluator;
import qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Atom.IBS;

export namespace Atoml
{
namespace Gaussian
{

class Orbital_IBS
    : public Gaussian_IBS
    , public ::Gaussian::IrrepBasisSet
    , public Atom::Orbital_HF_IBS <double>
    , public Atom::Orbital_IBS    <double>
    , public Atom::Orbital_DFT_IBS<double>
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
};

class Fit_IBS 
    : public Gaussian_IBS
    , public ::Gaussian::IrrepBasisSet
    , public Atom::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L);
    Fit_IBS(const DB_cache<double>* db,const ds_t& exponents, size_t L);
};

class BasisSet 
    : public Gaussian_BS 
    , public ::Gaussian::BS_Common
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
    BasisSet(const RVec& exponents, size_t Lmax); 
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
};

} //namespace
} //namespace


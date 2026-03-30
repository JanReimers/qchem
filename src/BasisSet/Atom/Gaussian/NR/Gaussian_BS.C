// File: BasisSet/Atom/Gaussian/NR/Gaussian_BS.C Gaussian Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.Gaussian.NR.BS;
import BasisSet.Atom.Gaussian.NR.IBS_EValuator;
import BasisSet.Atom.Gaussian_BS;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;

export namespace AtomBS
{
namespace Gaussian
{

class Orbital_IBS
    : public Gaussian_IBS
    , public AtomBS::IrrepBasisSet
    , public AtomBS::Orbital_HF_IBS <double>
    , public AtomBS::Orbital_IBS    <double>
    , public AtomBS::Orbital_DFT_IBS<double>
{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const rvec_t& exponents, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,const rvec_t& exponents, size_t L, const std::vector<int>& ml);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
};

class Fit_IBS 
    : public Gaussian_IBS
    , public AtomBS::IrrepBasisSet
    , public AtomBS::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,const rvec_t& exponents, size_t L);
};

class BasisSet 
    : public Gaussian_BS 
    , public ::BS_Common
    , public AtomIE_BS_2E<double> //HF support
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax);
    BasisSet(const rvec_t& exponents, size_t Lmax); 
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
private:
    void Insert(Orbital_IBS*);
};

} //namespace
} //namespace


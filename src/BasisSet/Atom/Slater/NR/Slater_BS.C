// File: BasisSet/Atom/Slater/NR/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.Slater.NR.BS;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.BS_Evaluator;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace Slater
{
class Orbital_IBS
    : public Slater_IBS
    , public Atom::IrrepBasisSet
    , public Atom::Orbital_HF_IBS <double>
    , public Atom::Orbital_IBS    <double>
    , public Atom::Orbital_DFT_IBS<double>
{
    using ds_t=std::valarray<double>;
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const ds_t& exponents, size_t L);
    Orbital_IBS(const DB_BS_2E<double>* db,const ds_t& exponents, size_t L, const std::vector<int>& ml);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
};

class Fit_IBS 
: public Slater_IBS
, public Atom::IrrepBasisSet
, public Atom::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,const ds_t& exponents, size_t L);
};

class BasisSet 
    : public Slater_BS 
    , public ::BS_Common
    , public AtomIE_BS_2E<double> //HF support
{
    using ds_t=std::valarray<double>;
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax); 
    BasisSet(const ds_t& exponents, size_t Lmax); 
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
private:
    void Insert(Orbital_IBS*);
};

}} //namespace Slater::l


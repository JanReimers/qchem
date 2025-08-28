// File: BasisSet/Atom/l/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.Internal.l.SlaterBS;
import BasisSet.Atom.Slater_IBS;
import BasisSet.Atom.Slater_BS;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;

export namespace Atoml
{
namespace Slater
{
class Orbital_IBS
    : public virtual IBS_Evaluator
    , public Slater_IBS
    , public Atom::IrrepBasisSet
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
: public Slater_IBS
, public Atom::IrrepBasisSet
, public Atom::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,const Vector<double>& exponents, size_t L);
    Fit_IBS(const DB_cache<double>* db,const ds_t& exponents, size_t L);
};

class BasisSet 
    : public Slater_BS 
    , public ::BS_Common
    , public AtomIE_BS_2E<double> //HF support
{
public:
    BasisSet(size_t N, double minexp, double maxexp, size_t Lmax); 
    BasisSet(const RVec& exponents, size_t Lmax); 
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
private:
    void Insert(Orbital_IBS*);
};

}} //namespace Slater::l


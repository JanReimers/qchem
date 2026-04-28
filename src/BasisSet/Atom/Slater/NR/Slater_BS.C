// File: BasisSet/Atom/Slater/NR/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
export module qchem.BasisSet.Atom.Slater.NR.BS;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.BS_Evaluator;
import qchem.BasisSet.Atom.IBS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.DB_Cache;

export namespace AtomBS
{
namespace Slater
{
class Orbital_IBS
    : public Slater_IBS
    , public AtomBS::IrrepBasisSet
    , public AtomBS::Orbital_HF_IBS <double>
    , public AtomBS::Orbital_IBS    <double>
    , public AtomBS::Orbital_DFT_IBS<double>
{
public:
    Orbital_IBS(const DB_BS_HF<double>* db,const rvec_t& exponents, const Irrep_QNs::sym_t&);
    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
};

class Fit_IBS 
: public Slater_IBS
, public AtomBS::IrrepBasisSet
, public AtomBS::Fit_IBS
{
public:
    Fit_IBS(const DB_cache<double>* db,const rvec_t& exponents, size_t L);
};

class BasisSet 
    : public Slater_BS 
    , public ::BS_Common
    , public AtomIE_BS_HF<double> //HF support
{
public:
    BasisSet(const rvec_t& exponents, const ElectronConfiguration& ec);
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
private:
    void Insert(Orbital_IBS*);
};

}} //namespace Slater::l


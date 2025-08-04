// File: Atom/ml/Slater_BS.H  r^l exp(-a*r)*Y_lm type basis set.
module;
#include <vector>

export module qchem.BasisSet.Atom.Internal.ml.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.l.SlaterBS;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisSet.Internal.HeapDB;
import qchem.HF_IBS;
import qchem.DFT_IBS;
import qchem.BasisSet.Atom.Internal.ml.Angular;
import qchem.BasisSet.Internal.Common;
import qchem.Symmetry.ElectronConfiguration;

export namespace Atom_ml
{
namespace Slater
{
class Orbital_IBS
    : public virtual TOrbital_HF_IBS<double>
    , public virtual TOrbital_DFT_IBS<double>
    , public         ::Slater::IrrepBasisSet
    , public         Orbital_IBS_Common1<double>
    , public         Orbital_DFT_IBS_Common<double>
    , public         Orbital_HF_IBS_Common<double>
    , public         Atoml::Slater::Orbital_IE

{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};
class BasisSet 
    : public ::Slater::BS_Common
    , public IE_BS_2E_Angular
{
public:
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
    
};

}} //namespace Slater_m


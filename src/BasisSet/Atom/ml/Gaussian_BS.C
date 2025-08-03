// File: Atom/ml/Gaussian_BS.H  r^l exp(-ar^2)*Y_lm type basis set.
module;
#include <vector>
export module qchem.BasisSet.Atom.Internal.ml.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Atom.Internal.l.GaussianBS;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisSet.Internal.IBS_Common1;
import qchem.BasisSet.Internal.HeapDB;
import qchem.HF_IBS;
import qchem.Fit_IBS;
import qchem.Cluster;
import qchem.BasisSet.Atom.Internal.ml.Angular;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.Common;
import qchem.Symmetry.ElectronConfiguration;

export namespace Atom_ml
{
namespace Gaussian
{
class BasisFunction : public Atoml::Gaussian::BasisFunction
    {
    public:
        BasisFunction(double theExponent,int n, int l, int ml, double norm);
        
        virtual BasisFunction* Clone(        ) const;
    private:
        typedef Atoml::Gaussian::BasisFunction Base;
        int ml;
    };
class Orbital_IBS
: public virtual TOrbital_HF_IBS<double>
// , public virtual TOrbital_DFT_IBS<double>
, public         ::Gaussian::IrrepBasisSet
, public         Orbital_IBS_Common1<double>
// , public         Orbital_DFT_IBS_Common<double>
, public         Orbital_HF_IBS_Common<double>
, public         Atoml::Gaussian::Orbital_IE

{
public:
    Orbital_IBS(const DB_BS_2E<double>* db,const Vector<double>& exponents, size_t L, const std::vector<int>& ml);

    virtual ::Fit_IBS* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
private:
    void InsertBasisFunctions();
};
class BasisSet 
: public ::Gaussian::BS_Common
, public IE_BS_2E_Angular //Pick angular integrals.
{
public:
    BasisSet() {};
    BasisSet(size_t N, double minexp, double maxexp, const ElectronConfiguration& ec);
};

} //namespace 
} //namespace 


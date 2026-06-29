// File: ChargeDensity/FourierSeedCD.C  Plane-wave (G-space) superposition-of-atomic-densities seed.
//
// The reciprocal-space twin of CompositeFittedCD: the SAD seed for the plane-wave (dcmplx) path.  Its native
// representation is rho-tilde(G) = Sum_atoms rho_atom(|G|) e^{-iG.R} -- a per-species radial form factor
// (the atomic VALENCE density's Fourier transform) assembled with structure factors, exactly the assembly
// the local pseudopotential already uses.  It is a tChargeDensity<dcmplx> (NOT a tDM_CD: a sum of atomic
// densities has no density matrix); the PW DFT terms (PW_Hartree, PW_XC) consume it through the
// FourierDensity face (GetFourierDensity), so no matrix is ever needed.
module;
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
export module qchem.ChargeDensity.FourierSeedCD;
export import qchem.ChargeDensity;                 // tChargeDensity<dcmplx>
export import qchem.ChargeDensity.FourierDensity;  // FourierDensity, FourierMap
import qchem.ChargeDensity.AtomicDensity;          // RadialDensity, RecentredAtomicDensity, GetAtomicDensity
import qchem.BasisSet.Band_FT_IBS;                  // Band_FT_IBS (the G-space assembler)
import qchem.Structure;                             // Structure, Atom
import qchem.ScalarFunction;                        // ScalarFunction<double>

export namespace qchem::ChargeDensity
{

class FourierSeedCD
    : public virtual tChargeDensity<dcmplx>
    , public virtual FourierDensity
{
public:
    //! Build the seed for plane-wave block \a basis and structure \a st: read each element's pseudo-valence
    //! radial density (\a functional, from atomic_valence_densities.json) and prepare the form-factor sum.
    FourierSeedCD(const BasisSet::Band_FT_IBS* basis, const Structure* st, const std::string& functional="LDA");

    // FourierDensity -- the native G-space representation the PW Hartree/XC terms consume.
    virtual FourierMap GetFourierDensity() const;

    // ScalarFunction<double> -- real-space rho(r) = Sum_atoms rho_atom(|r-R|) (not used by the PW Fock build,
    // which goes through GetFourierDensity, but provided so the type is whole).
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

    // tChargeDensity<dcmplx>
    virtual double GetTotalCharge() const {return itsScale*itsCharge;}
    virtual size_t Version       () const {return itsVersion;}
    virtual void   ReScale(double factor);

private:
    const BasisSet::Band_FT_IBS* itsBasis;        //!< the plane-wave block (G-space assembler); not owned
    const Structure*             itsStructure;    //!< atom Z + positions; not owned
    std::map<int,std::shared_ptr<const RadialDensity>> itsRadByZ;   //!< per-element valence radial density
    std::vector<RecentredAtomicDensity>                itsRecentred; //!< per-atom rho_atom(|r-R|) for op(r)
    double itsCharge;     //!< total valence electron count (Sum_atoms N_val), pre-scale
    double itsScale=1.0;  //!< uniform scale applied by ReScale
    size_t itsVersion;    //!< transient freshness serial
};

} //namespace

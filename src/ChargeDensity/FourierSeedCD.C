// File: ChargeDensity/FourierSeedCD.C  Plane-wave (G-space) superposition-of-atomic-densities seed.
//
// The reciprocal-space twin of NumericCD: the SAD seed for the plane-wave (dcmplx) path.  Its native
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
export import qchem.ChargeDensity.FourierDensity;  // FourierDensity, ΔG_Map
import qchem.ChargeDensity.AtomicDensity;          // RadialDensity, RecentredAtomicDensity, GetAtomicDensity
import qchem.BasisSet.Band_FT_IBS;                  // Band_FT_IBS (the G-space structure-factor assembler)
import qchem.BasisSet.Fit_IBS;                      // cFIT_CD_ABS (GetRepulsion3C's fit-basis arg)
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
    //! \a ionicScaleByZ is the per-species IonicSAD multiplier (empty => neutral SAD, all 1.0): species \c Z
    //! gets its valence density scaled by \c (N_val - q_Z)/N_val so the seed carries the formal charge \c q_Z
    //! (Na+ -> 0, F- -> 8/7); the total still integrates to the cell's electron count (sum of charges = 0).
    FourierSeedCD(const BasisSet::Band_FT_IBS* basis, const Structure* st, const std::string& functional="LDA",
                  const std::map<size_t,double>& ionicScaleByZ = {});

    // FourierDensity -- the native G-space representation the PW Hartree/XC terms consume.  The seed carries
    // no D, so its rho-tilde IS the structure-factor density (fit basis ignored); V_H applies the diagonal
    // kernel to it (bit-identical to the old Repulsion(rho-tilde) path).
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual ΔG_Map GetRepulsion3C   (const BasisSet::cFIT_CD_ABS& c) const;

    // ScalarFunction<double> -- real-space rho(r) = Sum_atoms rho_atom(|r-R|) (not used by the PW Fock build,
    // which goes through GetFourierDensity, but provided so the type is whole).
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

    // tChargeDensity<dcmplx>
    virtual double GetTotalCharge() const {return itsScale*itsCharge;}
    virtual size_t Version       () const {return itsVersion;}
    virtual void   ReScale(double factor);

private:
    ΔG_Map StructureFactorDensity() const;        //!< the seed's rho-tilde = per-species structure-factor sum
    const BasisSet::Band_FT_IBS* itsBasis;        //!< the plane-wave block (structure-factor assembler); not owned
    const Structure*             itsStructure;    //!< atom Z + positions; not owned
    std::map<size_t,std::shared_ptr<const RadialDensity>> itsRadByZ; //!< per-element (Z) valence radial density
    std::map<size_t,double>                            itsScaleByZ; //!< per-element IonicSAD multiplier (def 1.0)
    std::vector<RecentredAtomicDensity>                itsRecentred; //!< per-atom rho_atom(|r-R|) for op(r)
    double itsCharge;     //!< total (post-ionic-scale) valence electron count, pre-ReScale
    double itsScale=1.0;  //!< uniform scale applied by ReScale
    size_t itsVersion;    //!< transient freshness serial
};

} //namespace

// File: ChargeDensity/SeedCD.C  Plane-wave (G-space) superposition-of-atomic-densities seed.
//
// The reciprocal-space twin of NumericCD: the SAD seed for the plane-wave (dcmplx) path.  rho-tilde(G) =
// Sum_atoms F(Z,|G|) e^{-iG.R} -- the ANALYTIC structure-factor assembly of the per-species radial form
// factor F (an atomic VALENCE density's 1-D Fourier transform).  It builds this through its OWN density-fit
// basis (a cFIT_CD_ABS, which is a G_FieldEvaluator) -- NOT the orbital basis: the seed depends only on the
// fit basis it is handed.  (Grid-sampling rho(r)+FFT would ALIAS the peaked atomic density; the analytic form
// factor does not.)  It is a tChargeDensity<dcmplx> (NOT a tDM_CD: a sum of atomic densities has no density
// matrix); the PW DFT terms (PW_Hartree, PW_XC) consume it through the FourierDensity face.
module;
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
export module qchem.ChargeDensity.SeedCD;
export import qchem.ChargeDensity;                 // tChargeDensity<dcmplx>
export import qchem.ChargeDensity.FourierDensity;  // FourierDensity, ΔG_Map
import qchem.ChargeDensity.AtomicDensity;          // RadialDensity, RecentredAtomicDensity, GetAtomicDensity
import qchem.BasisSet.Fit_IBS;                      // cFIT_CD_ABS (the density-fit basis it builds rho-tilde through)
import qchem.Structure;                             // Structure, Atom
import qchem.ReciprocalLattice;                     // ReciprocalLattice (the seed's own Poisson metric B)
import qchem.ScalarFunction;                        // ScalarFunction<double>

export namespace qchem::ChargeDensity
{

class SeedCD
    : public virtual tChargeDensity<dcmplx>
    , public virtual FourierDensity
{
public:
    //! Build the seed for density-fit basis \a fitBasis (from the orbital basis's CreateCDFitBasisSet) and
    //! structure \a st: read each element's pseudo-valence radial density (\a functional, from
    //! atomic_valence_densities.json) and prepare the form-factor sum.  \a ionicNvalByZ is the per-species
    //! IonicSAD TARGET valence-electron count (empty => neutral SAD): species \c Z should carry \c Nval-q_Z
    //! electrons (Na+ -> 0, F- -> 8).  For each species we prefer the library's CHARGE-STATE density with that
    //! \c Nelec (a DIFFUSE F- anion -- what makes the ionic seed actually converge), scale 1; if the library
    //! lacks it we fall back to the neutral density amplitude-scaled by \c target/Nval (the old compact seed);
    //! a target of 0 (a stripped cation) is the zero density.  The total integrates to the cell electron count.
    SeedCD(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitBasis, const Structure* st,
                  const std::string& functional="LDA", const std::map<size_t,int>& ionicNvalByZ = {});

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
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> itsFitBasis; //!< density-fit basis (its grid engine builds rho-tilde); owned
    const Structure*             itsStructure;    //!< atom Z + positions; not owned
    ReciprocalLattice            itsRecip;        //!< the cell's reciprocal lattice (B): the seed's Poisson metric
    std::map<size_t,std::shared_ptr<const RadialDensity>> itsRadByZ; //!< per-element (Z) valence radial density
    std::map<size_t,double>                            itsScaleByZ; //!< per-element IonicSAD multiplier (def 1.0)
    std::vector<RecentredAtomicDensity>                itsRecentred; //!< per-atom rho_atom(|r-R|) for op(r)
    double itsCharge;     //!< total (post-ionic-scale) valence electron count, pre-ReScale
    double itsScale=1.0;  //!< uniform scale applied by ReScale
    size_t itsVersion;    //!< transient freshness serial
};

} //namespace

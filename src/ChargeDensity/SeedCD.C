// File: ChargeDensity/SeedCD.C  Plane-wave (G-space) superposition-of-atomic-densities seed.
//
// The reciprocal-space twin of NumericCD: the SAD seed for the plane-wave (dcmplx) path.  rho-tilde(G) =
// Sum_atoms F(Z,|G|) e^{-iG.R} -- the ANALYTIC structure-factor assembly of the per-species radial form
// factor F (an atomic VALENCE density's 1-D Fourier transform).  It builds this through its OWN density-fit
// basis (a cFIT_CD_ABS, which is a G_FieldEvaluator) -- NOT the orbital basis: the seed depends only on the
// fit basis it is handed.  (Grid-sampling rho(r)+FFT would ALIAS the peaked atomic density; the analytic form
// factor does not.)  It is a tChargeDensity<dcmplx> (NOT a tDM_CD: a sum of atomic densities has no density
// matrix); the PW DFT terms (PW_Hartree, the XC quadrature) consume it through the FourierDensity face.
//
// SPIN (doc/SCFSeedingPlan.md §10): a SeedCD can also be ONE CHANNEL of the spin-polarized seed -- the
// \c channel ctor argument selects Up/Down, reading the library's spin-resolved (majority, minority) pair
// per species (up-majority storage) and honouring each atom's assembly-time \c itsSpinFlip bit (the -m
// sublattice swaps the pair -- the AFM staggering).  PolarizedSeedCD below composes two channel SeedCDs
// into the two-channel seed the polarized SCF consumes.
module;
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
export module qchem.ChargeDensity.SeedCD;
export import qchem.ChargeDensity;                 // tChargeDensity<dcmplx>, cSpinResolved_CD, Spin
export import qchem.ChargeDensity.FourierDensity;  // FourierDensity, ΔG_Map
import qchem.ChargeDensity.AtomicDensity;          // RadialDensity, RecentredAtomicDensity, GetAtomicDensity
import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the density-fit basis it builds rho-tilde through)
import qchem.Structure;                             // Structure, Atom
import qchem.UnitCell;                              // UnitCell (the flip-group sub-cells of a channel seed)
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
    //!
    //! \a channel (default \c Spin::None = the historical spin-agnostic total) makes this ONE CHANNEL of the
    //! polarized seed: a species with a library spin pair AT THE TARGET \c Nelec contributes its majority
    //! (\c Up) / minority (\c Down) component -- swapped on atoms with \c itsSpinFlip set (the AFM -m
    //! sublattice) -- and a species without a pair contributes the spin-agnostic half density (rho/2, both
    //! channels alike; the flip is then a no-op).  \c None ignores flips entirely (the total is flip-blind).
    SeedCD(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitBasis, const Structure* st,
                  const std::string& functional="LDA", const std::map<size_t,int>& ionicNvalByZ = {},
                  const Spin& channel = Spin::None);

    // FourierDensity -- the native G-space representation the PW Hartree/XC terms consume.  The seed carries
    // no D, so its rho-tilde IS the structure-factor density (fit basis ignored); V_H applies the diagonal
    // kernel to it (bit-identical to the old Repulsion(rho-tilde) path).
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual ΔG_Map GetRepulsion3C   (const BasisSet::cFIT_CD_ABS& c) const;

    // ScalarFunction<double> -- real-space rho(r) = Sum_atoms rho_atom(|r-R|) (not used by the PW Fock build,
    // which goes through GetFourierDensity, but provided so the type is whole).
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;
    //! Batch \f$\rho(r_g)\f$ over a whole mesh -- the \c ScalarFunction default loop, THREADED (opt-in,
    //! \c GPW_OMP_THREADS).  This is how the atom-centred XC route samples the SEED (the one iteration
    //! with no density matrix to GEMM), and on a magnetic cell it was the largest single bucket of the
    //! whole run: 28 s of image-summed radial lookups over 48k mesh points x 2 channels.  Every point is
    //! independent and \c op() is pure, so the raster is bit-identical at any thread count.
    virtual rvec_t  operator()(const rvec3vec_t&) const;

    // tChargeDensity<dcmplx>
    virtual double GetTotalCharge() const {return itsScale*itsCharge;}
    virtual size_t Version       () const {return itsVersion;}
    virtual void   ReScale(double factor);

private:
    ΔG_Map StructureFactorDensity() const;        //!< the seed's rho-tilde = per-species structure-factor sum
    //! One flip-group's structure-factor sum: \a group's atoms against \a radByZ's per-species radial.
    ΔG_Map GroupDensity(const Structure* group, const std::map<size_t,std::shared_ptr<const RadialDensity>>& radByZ) const;
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> itsFitBasis; //!< density-fit basis (its grid engine builds rho-tilde); owned
    const Structure*             itsStructure;    //!< atom Z + positions (+ spin-flip bits); not owned
    ReciprocalLattice            itsRecip;        //!< the cell's reciprocal lattice (B): the seed's Poisson metric
    const UnitCell*              itsCell;         //!< == itsStructure (checked at construction): the direct cell, for the periodic real-space image sums
    Spin                         itsChannel;      //!< None = spin-agnostic total; Up/Down = one polarized channel
    std::map<size_t,std::shared_ptr<const RadialDensity>> itsRadByZ;     //!< per-element channel radial, UNFLIPPED sites
    std::map<size_t,std::shared_ptr<const RadialDensity>> itsRadFlipByZ; //!< per-element channel radial, FLIPPED sites (== itsRadByZ unless a spin pair)
    std::map<size_t,double>                            itsScaleByZ; //!< per-element multiplier (IonicSAD fallback x the pairless rho/2)
    //! Flip-partitioned sub-cells for the G assembly (built only when a channel seed has flipped atoms --
    //! the basis's MakeFourierDensity is species-keyed, so each flip group gets its own call).
    std::shared_ptr<UnitCell>                          itsGroupA, itsGroupB;
    std::vector<RecentredAtomicDensity>                itsRecentred; //!< per-atom rho_atom(|r-R|) for op(r), flip-aware
    std::vector<double>                                itsScalePerAtom; //!< per-atom multiplier, parallel to itsRecentred (keeps op(r) off the structure)
    double itsCharge;     //!< total (post-ionic-scale) valence electron count, pre-ReScale
    double itsScale=1.0;  //!< uniform scale applied by ReScale
    size_t itsVersion;    //!< transient freshness serial
};

//! The SPIN-POLARIZED plane-wave SAD seed (doc/SCFSeedingPlan.md §10): two channel SeedCDs (majority
//! tables + per-atom flip bits = the collinear magnetic configuration, e.g. the MnO AFM-II staggering).
//! The matrix-free polarized sibling of cPolarized_CD: the spin-native XC engine reads the channels
//! through the cSpinResolved_CD face; the total-density consumers (PW_Hartree) see the ↑+↓ sum through
//! the same FourierDensity face every seed has.  The polarized twin of NumericCD/SeedCD -- no density
//! matrix, so none of the tDM_CD contract appears anywhere.
class PolarizedSeedCD
    : public virtual tChargeDensity<dcmplx>
    , public virtual FourierDensity
    , public virtual cSpinResolved_CD
{
public:
    //! Same arguments as SeedCD; builds the Up + Down channel pair (flip bits read from \a st's atoms).
    PolarizedSeedCD(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitBasis, const Structure* st,
                    const std::string& functional="LDA", const std::map<size_t,int>& ionicNvalByZ = {});

    // cSpinResolved_CD -- the spin-native XC engine's channel access.
    virtual const cChargeDensity* GetChannel(const Spin&) const;

    // FourierDensity -- the ↑+↓ TOTAL (what Hartree consumes; spin never enters V_H).
    virtual ΔG_Map GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const;
    virtual ΔG_Map GetRepulsion3C   (const BasisSet::cFIT_CD_ABS& c) const;

    // ScalarFunction<double> -- the total rho(r).
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;

    // tChargeDensity<dcmplx>
    virtual double GetTotalCharge() const;
    double         GetTotalSpin  () const;   //!< <up>-<down>: ~0 for an AFM staggering, 2S for FM/Hund
    virtual size_t Version       () const {return itsVersion;}
    virtual void   ReScale(double factor);

private:
    SeedCD itsUp, itsDn;   //!< the two channel seeds (share tables via the library; cheap value members)
    size_t itsVersion;     //!< transient freshness serial (one serial for the pair)
};

} //namespace

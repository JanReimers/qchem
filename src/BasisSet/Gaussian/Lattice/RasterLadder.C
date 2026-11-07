// File: BasisSet/Gaussian/Lattice/RasterLadder.C  THE GPW RASTER LADDER: the density grid and its REL_CUTOFF multigrid levels, as ONE
// object (V-BS-RASTER, user + controller 2026-11-02).
//
// A basis set has TWO jobs (user): serve its irrep blocks, and create the CD / XC fit bases.  The raster ladder is neither, yet every
// k-block built its own copy (GPW_IBS IS-A GPW_Evaluator, whose ctor built the density grid and EnsureLevels the multigrid) -- from
// k-INDEPENDENT inputs: the cell, the molecular basis's exponent range, densityEcut / cutoffFactor / ladderFactor, the raster-fields
// routing policy and the GPWTolerances.  So a k222 run held eight identical ladders.  This is the ladder built ONCE per geometry point
// (the cell moves in a SolidSession) and shared, by \c MakeRasterLadder.
//
// THE BASIS INPUT (step 2, 2026-11-07) is the basis's own answer to "how fine must a grid over your density reach", \c MakeResolution
// at \c kResolutionEps: the ladder reads the sharpest and the most diffuse pair product off it (the basis owns the factor 2), and asks
// no exponent of any block.  The pair rule's stiffness is the run's \c GPWTolerances::relCutoffSafety (a named calibration), the
// same value the basis's own pair->level assignment is bound with.
// MULTIGRID-READY: the pair -> level cap (\c PairCap) is the ladder's, so a per-primitive level scheme replaces only its internals.
module;
#include <cstddef>
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
export module qchem.BasisSet.Gaussian.Lattice.RasterLadder;
export import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // the typed tolerances the ladder is sized by
export import qchem.BasisSet.PlaneWave.Evaluators;             // PlaneWave::PW_Grid_Evaluator (a level) + RasterFields
export import qchem.BasisSet.GridResolution;                     // Resolution: the basis's answer the ladder is sized from
import qchem.BasisSet;                                           // Real_BS (the molecular basis a ladder is sized from)
import qchem.BasisSet.Orbital_PP_IBS;                            // SpeciesRadialField(_Gaussian): the local PP's sharpness
import qchem.Symmetry.Lattice_3D.SpaceGroup;                     // DirectOp {W|τ} (the grids report's {G}-star column)
import qchem.Structure;
import qchem.UnitCell;
import qchem.ReciprocalLattice;
import qchem.Types;

export namespace qchem::BasisSet::Gaussian
{

//! The k-independent statement a ladder is built from (beside the cell and the basis): what a caller that hands a ladder in must
//! have meant, so a block can refuse a ladder built for other settings instead of silently running on it.
struct RasterLadderInputs
{
    double        densityEcut  = -1.0;
    double        cutoffFactor = 2.0;
    double        ladderFactor = 4.0;
    PlaneWave::RasterFields rasterFields = PlaneWave::RasterFields::HartreeXC;
    GPWTolerances tol = {};
    bool operator==(const RasterLadderInputs&) const = default;
};

//! One ladder over a fine grid: level [0] is that grid (the resolution reference), coarser levels follow, the top completion rung LAST.
struct RasterLevels
{
    std::vector<std::shared_ptr<const PlaneWave::PW_Grid_Evaluator>> levels;
    std::vector<ivec3_t> N;          //!< each level's FFT grid divisions
    std::vector<double>  ecut;       //!< each level's cutoff (reference first)
    std::size_t          nBase = 0;  //!< levels before the top rung (the local-PP sub-ladder)
};

//! \brief The LOCAL pseudopotential's sharpness \f$\beta\f$ over the atoms of \a cl: V_loc long's effective Gaussian exponent, the
//! SHARPEST species' \f$1/(2 r_{loc}^2)\f$ (V_long sums over atoms) -- the field side of V_long's harmonic routing rule.  0 when \a loc
//! has no closed Gaussian form (V_long then takes the kappa-sweep path and has no private ladder).  ONE home: \c MakeLocalPPLong and the
//! facade that sizes the ladder both call it, so the ladder's V_long audit and the term's ladder agree by construction.
double LocalPPSharpness(const Structure& cl, const SpeciesRadialField& loc);

//! \brief The GPW density grid and its multigrid ladder, built once from k-independent inputs and shared by every block of a point.
class GpwRasterLadder
{
public:
    //! \param cell the direct lattice (its reciprocal cell spans every level's {G}).
    //! \param resolution the molecular basis's \c MakeResolution at \a resolutionEps: its density's finest and coarsest cutoffs, from
    //!        which the ladder reads the sharpest and the most diffuse pair product (\f$p=\tfrac12G^2/\ln(1/\epsilon)\f$).
    //! \param densityEcut <0 AUTO (floor \a cutoffFactor * alpha_max), 0 no DFT tier (no grid), >0 explicit Ha (warned if under the floor).
    //! \param ppSharpness the local PP's \c LocalPPSharpness (0: none known -- a bare block's ladder; V_long then files no audit row
    //!        here).  STEP 1 states it as this one number; step 2 replaces it with the PP set's own resolution answer, combined with
    //!        the basis's (the PP sharper than the orbital basis is exactly what it must size).
    GpwRasterLadder(const UnitCell& cell, const BasisSet::Resolution& resolution, double resolutionEps, double densityEcut, double cutoffFactor,
                    double ladderFactor, PlaneWave::RasterFields rasterFields, const GPWTolerances& tol, double ppSharpness=0.0);

    //! The tolerance a ladder asks the basis's resolution at.  Any eps in (0,1) names the same pair exponents (the Gaussian answer is
    //! linear in ln(1/eps)); this one is the collocation floor's, so the round trip sits at the ulp.
    static constexpr double kResolutionEps = 1.0e-10;

    //! The density / collocation grid (level 0 of \c Levels()); null when the DFT tier is off (densityEcut 0).
    const std::shared_ptr<const PlaneWave::PW_Grid_Evaluator>& DensityGrid() const {return itsGrid;}
    //! The ladder over the density grid (built with the ladder, never lazily).  \pre \c DensityGrid() is set.
    const RasterLevels& Levels() const;
    //! The ladder over a REQUESTED fine grid (a fit basis's): the same rule, a different reference.
    RasterLevels MakeLevels(std::shared_ptr<const PlaneWave::PW_Grid_Evaluator> grid) const;
    //! The tightest pair exponent a level at \a ecut may carry (\a finest: every pair), \a ecutRef the reference level's cutoff.
    double PairCap(double ecut, double ecutRef, bool finest) const;
    //! The HartreeOnly routing floor (a fraction of alpha_max), or <0 when the run routes every pair at the XC floor (HartreeXC).
    double RelFieldSharp() const {return itsRelFieldSharp;}
    //! \brief V_loc long's PRIVATE ladder (\c GPW_Evaluator::MakeLocalPPLong): the density ladder plus ONE top level at the sharpest
    //! pair's harmonic requirement \f$2 L\,p\beta/(p+\beta)\f$ (\f$p=2\alpha_{\max}\f$, \a lnE \f$=L=\ln(1/\epsilon)\f$) when no level
    //! reaches it, every level re-sampled at V_long's own raster demand when it exceeds the density's.  \pre \c DensityGrid() is set.
    RasterLevels MakeVlocLongLevels(double beta, double lnE) const;

    //! The local PP sharpness the ladder was sized with (0: none) -- a block's V_long checks its own field against it.
    double PPSharpness() const {return itsPPSharpness;}

    // --- THE LADDER REPORTS ITS OWN GRIDS (V-BS-RASTER): every answer below is a pure function of the ladder, so the facade ASKS it
    //     (never a block, never the report) and no row depends on which term happened to run.
    //! The alias exponents the rasters deliver (V_long's from its private ladder at \c PPSharpness; 0 without one).  None (all 0)
    //! when the DFT tier is off.
    RasterExponents DeliveredRasterExponents() const;
    //! The grid-audit verdicts: the density raster's, then V_long's (with a \c PPSharpness).  Empty when the DFT tier is off.
    std::vector<RasterVerdict> DeliveredRasterVerdicts() const;
    //! \brief The `grids` facts, written CURSOR-RELATIVE into the caller's open scope (the facade opens "grids"): the density-floor
    //! settings, the ladder rows (with the {G}-star count under \a ops -- interim, V-BS-SYMMETRY moves it), the local-PP and V_long
    //! settings, and the audit rows (\c auditDensity, \c auditVlocLong), each with its UNDER line on the console when it fails.
    void ReportGrids(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops) const;
    //! The same facts as console lines, for a run with no report open (\a nFunctions: the basis size the header line names).
    void PrintGrids(std::ostream& os, size_t nFunctions, const std::vector<Symmetry::Lattice_3D::DirectOp>& ops) const;

    const RasterLadderInputs& Inputs() const {return itsInputs;}   //!< what the ladder was built from (a handed-in ladder is checked)
    double CutoffFactor() const {return itsInputs.cutoffFactor;}   //!< C in the auto density floor C*alpha_max
    double LadderFactor() const {return itsInputs.ladderFactor;}   //!< the coarse-level progression factor

private:
    std::shared_ptr<const PlaneWave::PW_Grid_Evaluator> MakeLevel(const ReciprocalLattice& recip, double ecut, double pCap, double x) const;

    RasterLadderInputs itsInputs;
    UnitCell      itsCell;          //!< the direct cell (the audits' raster Nyquist radius)
    double        itsPPSharpness=0.0;
    //! The sharpest and the most diffuse FUNCTION exponents, read off the basis's resolution: half its sharpest / most diffuse pair
    //! product.  The ladder's arithmetic (the floor C*alpha_max, the pair caps, the field floors) is written in them.
    double        itsAlphaMax, itsAlphaMin;
    double        itsRelCutoffSafety;   //!< the run's GPWTolerances::relCutoffSafety (EMPIRICAL, record at kRelCutoffSafety)
    double        itsLadderFactor;
    double        itsRelFieldSharp=-1.0;
    GPWTolerances itsTol;
    std::shared_ptr<const PlaneWave::PW_Grid_Evaluator> itsGrid;
    RasterLevels  itsLevels;
};

//! \brief Build a ladder over \a cell from \a mol's one Gaussian orbital block (its \c MakeResolution).  The one place a ladder is made
//! from a molecular basis: \c BasisSet::Lattice::MakeRasterLadder (per geometry point) and a bare block's convenience ctor both call it.
//! Throws if \a mol has more than one orbital block, or none that is Gaussian.
std::shared_ptr<const GpwRasterLadder> MakeRasterLadder(const UnitCell& cell, const BasisSet::Real_BS& mol, const RasterLadderInputs& in,
                                                        double ppSharpness=0.0);

//! One diagnostic line for any stored grid: FFT divisions, Ecut, nG and the spectral extent |G|min / |G|max (also printed by the fit
//! factories for the grid they wrap).
std::ostream& PrintGrid(std::ostream& os, const std::string& tag, const PlaneWave::PW_Grid_Evaluator& g);

} // namespace qchem::BasisSet::Gaussian

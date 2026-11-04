// File: BasisSet/Gaussian/Lattice/RasterLadder.C  THE GPW RASTER LADDER: the density grid and its REL_CUTOFF multigrid levels, as ONE
// object (V-BS-RASTER, user + controller 2026-11-02).
//
// A basis set has TWO jobs (user): serve its irrep blocks, and create the CD / XC fit bases.  The raster ladder is neither, yet every
// k-block built its own copy (GPW_IBS IS-A GPW_Evaluator, whose ctor built the density grid and EnsureLevels the multigrid) -- from
// k-INDEPENDENT inputs: the cell, the molecular basis's exponent range, densityEcut / cutoffFactor / ladderFactor, the raster-fields
// routing policy and the GPWTolerances.  So a k222 run held eight identical ladders.  This is the ladder built ONCE per geometry point
// (the cell moves in a SolidSession) and shared, by \c MakeRasterLadder.
//
// STEP 1 (bit-identical): the ladder's basis input is still the EXPONENT RANGE (alpha_max, alpha_min, RelCutoffSafety) through
// GaussianSharpness, exactly the arithmetic the evaluator does.  Step 2 switches it to \c MakeResolution(eps), judged by a sweep.
// MULTIGRID-READY: the pair -> level cap (\c PairCap) is the ladder's, so a per-primitive level scheme replaces only its internals.
module;
#include <cstddef>
#include <memory>
#include <vector>
export module qchem.BasisSet.Gaussian.Lattice.RasterLadder;
export import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // the typed tolerances the ladder is sized by
export import qchem.BasisSet.PlaneWave.Evaluators;             // PlaneWave::PW_Grid_Evaluator (a level) + RasterFields
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;              // GaussianSharpness: the exponent range (step 1's basis input)
import qchem.BasisSet;                                           // Real_BS (the molecular basis a ladder is sized from)
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

//! \brief The GPW density grid and its multigrid ladder, built once from k-independent inputs and shared by every block of a point.
class GpwRasterLadder
{
public:
    //! \param cell the direct lattice (its reciprocal cell spans every level's {G}).
    //! \param basis the molecular Gaussian basis's exponent range (step 1; \c MakeResolution in step 2).
    //! \param densityEcut <0 AUTO (floor \a cutoffFactor * alpha_max), 0 no DFT tier (no grid), >0 explicit Ha (warned if under the floor).
    GpwRasterLadder(const UnitCell& cell, const GaussianSharpness& basis, double densityEcut, double cutoffFactor, double ladderFactor,
                    PlaneWave::RasterFields rasterFields, const GPWTolerances& tol);

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

    const RasterLadderInputs& Inputs() const {return itsInputs;}   //!< what the ladder was built from (a handed-in ladder is checked)
    double CutoffFactor() const {return itsInputs.cutoffFactor;}   //!< C in the auto density floor C*alpha_max
    double LadderFactor() const {return itsInputs.ladderFactor;}   //!< the coarse-level progression factor

private:
    std::shared_ptr<const PlaneWave::PW_Grid_Evaluator> MakeLevel(const ReciprocalLattice& recip, double ecut, double pCap, double x) const;

    RasterLadderInputs itsInputs;
    double        itsAlphaMax, itsAlphaMin, itsRelCutoffSafety;
    double        itsLadderFactor;
    double        itsRelFieldSharp=-1.0;
    GPWTolerances itsTol;
    std::shared_ptr<const PlaneWave::PW_Grid_Evaluator> itsGrid;
    RasterLevels  itsLevels;
};

//! \brief Build a ladder over \a cell from \a mol's one Gaussian orbital block (its exponent range).  The one place a ladder is made
//! from a molecular basis: \c BasisSet::Lattice::MakeRasterLadder (per geometry point) and a bare block's convenience ctor both call it.
//! Throws if \a mol has more than one orbital block, or none that is Gaussian.
std::shared_ptr<const GpwRasterLadder> MakeRasterLadder(const UnitCell& cell, const BasisSet::Real_BS& mol, const RasterLadderInputs& in);

} // namespace qchem::BasisSet::Gaussian

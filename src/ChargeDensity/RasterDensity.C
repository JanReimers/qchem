// File: ChargeDensity/RasterDensity.C  A periodic density's raw ρ_DM(r) on a fit basis's raster: the INTERIM face.
//
// V-ONE-OF-N 3c (doc/Records/VOneOfN_D2Design.md §3) retired the FourierDensity face this file used to hold: ρ̃(G) is now asked
// through tChargeDensity::Project onto a BasisSet::GSpaceDensityProjector, the Hartree projection through the Hartree fitter, and
// the star-average happens inside the k-composite's Project.  What is left is the one capability the design note keeps for
// sub-step 4 (R1.0n/o): the density's RAW ρ_DM(r) on the raster, the non-negative grid the XC shadow and the ρ̃ mixers read.
// It is still a cross-cast face with an empty default ("no raw path": a plane-wave basis's collocation has none), the shape it
// had; sub-step 4 makes it a Projections realization and retires this file.
module;
export module qchem.ChargeDensity.RasterDensity;
import qchem.BasisSet.Orbital_DFT_IBS;   // cFIT_SF_ABS (the raster's fit basis)
import qchem.Types;                      // rvec_t

export namespace qchem::ChargeDensity
{

//! \brief A periodic density that may answer its raw \f$\rho_{DM}(\vec r)\f$ on a fit basis's raster.  Data-free.
class RasterDensity
{
public:
    virtual ~RasterDensity() = default;
    //! \brief Raw \f$\rho_{DM}=\sum_k w_k\rho_k(\vec r)\f$ on \a c's raster, from the basis's collocation-native forward (GPW's
    //! \c applyRaw); EMPTY when there is no raw path (plane waves), so a caller falls back to the ball route.  A composite answers
    //! empty unless every block answers (all or nothing), and star-averages in real space.
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS&) const {return rvec_t{};}
};

//! The empty stand-in a FINITE composite inherits instead (no raster), keeping its layout unchanged.
struct NoRasterDensity {};

} //namespace

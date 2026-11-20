// File: BasisSet/RasterDensityProjector.C  A density's RAW ρ_DM on a plane-wave fit basis's raster: the δ projector on the raster points.
//
// R1.0n/o (lane S 2026-11-19, controller-accepted design lane-s-V139-R10no-design.md).  The GPW XC feed and the ρ̃ mixers' raster
// shadow want \f$\rho_{DM}(r_g)\f$ on a plane-wave fit basis's uniform raster, collocated RAW (out-of-ball content included: pin 28).
// Those values are \f$\langle\delta_g|\rho\rangle\f$, the density's projections onto δ functions at the raster points -- so they are
// asked the one way a density answers anything, \c tChargeDensity::Project, onto THIS projector, and come back as \c Projections
// holding an \c rvec_t.  Until 2026-11-19 they were a second entry beside it (\c RasterDensity::GetRhoOnGrid, empty meaning "none"),
// which every client cross-cast to and tested for emptiness.
//   - a density with a matrix contracts its D into the block's collocation-native RAW forward (\c PairProjection::applyRaw);
//   - a density with no D (the mixers' presentation, the SAD seed) is asked for its field's raster values (\c RasterValuedField);
//   - a k-composite star-averages the sum once, IN REAL SPACE (a voxel permutation): the projector says how, the composite's group
//     averages -- the twin of the ρ̃ projector's G-space star-average.
// WHETHER a run has the raw route is decided ONCE, at setup, from the ORBITAL lineage (\c RasterPairCollocation): a Gaussian (GPW)
// block collocates its pair products on a raster, a plane-wave block does not.  A projector asked to forward a block that cannot
// THROWS: the setup that built it made a composition error.
module;
#include <memory>
#include <string>
export module qchem.BasisSet.RasterDensityProjector;
export import qchem.BasisSet.PairProjector;     // what the factory returns: the one projection face a density answers
export import qchem.BasisSet.Orbital_DFT_IBS;   // cFIT_SF_ABS (the raster's fit basis)
export import qchem.Types;                      // rvec_t

export namespace qchem::BasisSet
{

//! \brief "I collocate my orbital pair products RAW on a fit basis's raster": the orbital-lineage capability the setup asks before it
//! builds a raster projector (a GPW block has it; a plane-wave block, whose pair products ARE plane waves, does not).  Data-free: having
//! the face is the answer, and the projector's forward then reads the block's own raw collocation.
class RasterPairCollocation
{
public:
    virtual ~RasterPairCollocation() = default;
};

//! \brief A field with NO density matrix that carries its own values on a raster: the ρ̃ mixers' presentation (its raw-raster shadow)
//! and the SAD seed.  The raster projector asks it with the usual cross-cast between abstract faces.  Data-free.
class RasterValuedField
{
public:
    virtual ~RasterValuedField() = default;
    //! My values on \a raster's points, in raster order.  EMPTY when this field carries none: a seed has no raw values, and a mixer whose
    //! raw pipeline is off (or not yet activated) carries no shadow -- the mixer's own state, which it reads back as "off".
    virtual rvec_t RasterValues(const cFIT_SF_ABS& raster) const=0;
};

//! \brief The raw ρ_DM projector on the plane-wave fit basis \a sf's raster.  Clients hold the FACE and only ever ask
//! \c density.Project(*p).  THROWS when \a sf has no raster (\c G_RasterTransform) or cannot star-average one (\c G_RasterStarField).
std::unique_ptr<PairProjector<dcmplx>> MakeRasterDensityProjector(std::shared_ptr<const cFIT_SF_ABS> sf);

} //namespace

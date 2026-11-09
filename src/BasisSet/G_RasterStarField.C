// File: BasisSet/G_RasterStarField.C  A raster that knows how its voxels transform under a crystal op.
//
// The raster half of the density star-average (V-SYM-SEQUESTER-OPS, first piece; user rulings 2026-11-08/09).  The k-composite
// holds the GROUP and asks it to average its summed raster ρ_DM; the group iterates its own ops, and the raster supplies the
// image of its field under one op -- a voxel permutation plus the exact FFT glide shift.  It replaces
// G_RasterTransform::Symmetrize (removed from that tier-1 hub, user 2026-11-09), which averaged over ops the FIT BASIS stored:
// two instances of one basis then behaved differently, which an ID-keyed cache could not see.  The basis now stores no ops.
module;
#include <memory>
export module qchem.BasisSet.G_RasterStarField;
export import qchem.Symmetry.StarGroup;   // StarField
import qchem.Types;                       // rvec_t

export namespace qchem::BasisSet
{

//! \brief A raster engine's capability: wrap a field on my voxels (\c G_RasterTransform::RasterSize values, raster order) as
//! a \c Symmetry::StarField, whose image under \f$\{W|\tau\}\f$ is \f$f(W\cdot x+\tau)\f$.  Reached by a cross-cast from the
//! fit face that produced the raster (abstract to abstract).  Data-free.
class G_RasterStarField
{
public:
    virtual ~G_RasterStarField() = default;
    //! The star field over \a f, which it averages IN PLACE and must not outlive.
    virtual std::unique_ptr<Symmetry::StarField> MakeStarField(rvec_t& f) const=0;
};

} // namespace qchem::BasisSet

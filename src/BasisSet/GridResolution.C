// File: BasisSet/GridResolution.C  "How sharp is your sharpest function and how diffuse your most diffuse": the capability face a
// GRID BUILDER asks a basis (V-BS-RASTER, user 2026-11-02).
//
// A grid that samples a basis's density (a GPW raster ladder, a molecular radial mesh) must reach as far in G as the density's
// sharpest pair product and as far out in r as its most diffuse function.  Until now every builder asked for EXPONENTS
// (GaussianSharpness::MaxExponent / MinExponent) and applied its own factor of 2 for the pair product, so the question was Gaussian
// and the arithmetic was the caller's.  This face asks it in the grid's own language -- energies and lengths at a tolerance -- so a
// plane-wave or radial block can answer too, and the basis (not every caller) owns the step from its functions to its density.
//
// A TAIL AT A TOLERANCE, not an expectation value (doc/Pins.md pin 1: there is no cut in r space -- this says where a series may stop
// at eps, it imposes no radius on anything).
module;
#include <limits>
export module qchem.BasisSet.GridResolution;

export namespace qchem::BasisSet
{

//! What a grid builder needs from a basis at tolerance \a eps (a.u.: Ha, bohr).  \a eps is RELATIVE: each tail is measured against
//! its own PEAK amplitude (\f$G=0\f$ for a spectrum, \f$r=0\f$ for a function), never an absolute value.
//! \note A Gaussian block answers from its exponents, i.e. for s primitives.  An \f$l>0\f$ primitive carries an \f$r^l\f$ (and its
//! spectrum a \f$G^l\f$) prefactor, so its tail at the returned distance sits ABOVE eps by at most \f$(2eL/l)^{l/2}\f$,
//! \f$L=\ln(1/\epsilon)\f$ (the peak of \f$r^le^{-\alpha r^2}\f$ is at \f$r_0^2=l/2\alpha\f$) -- a few decades for d at
//! \f$\epsilon=10^{-10}\f$, i.e. ~10-20% more cutoff.  Folding it in needs each exponent's l, which the grid builders have not
//! asked for (UTGaussian_BS GridResolution.AnLShellsTailStaysWithinItsPrefactorMargin pins the bound).
struct Resolution
{
    //! \f$\tfrac12G^2\f$ beyond which the basis's DENSITY -- its pair products, so twice the sharpest exponent -- has Fourier
    //! amplitude below eps: the FINEST level a density grid needs.
    double cutoff       = 0.0;
    //! The same for the most diffuse pair product: the COARSEST level that still resolves any part of the density (a ladder's depth
    //! is the ratio of the two).
    double coarseCutoff = 0.0;
    //! The distance beyond which the most diffuse function is below eps (infinite for a delocalized block).
    double radius       = std::numeric_limits<double>::infinity();
};

//! \brief The face: a block that can say how fine and how far a grid over it must reach.
//!
//! Clients and what they do with the answer: the GPW raster ladder sizes its finest level from \c cutoff and its depth from
//! \c coarseCutoff / \c cutoff; the automatic density cutoff floors itself on \c cutoff; a molecular radial mesh may size its extent
//! from \c radius.  Reduced over a whole basis by \c tBasisSet::MakeResolution.
class GridResolution
{
public:
    virtual ~GridResolution() = default;
    //! \pre 0 < \a eps < 1.
    virtual Resolution MakeResolution(double eps) const = 0;
};

} // namespace qchem::BasisSet

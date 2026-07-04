//! \file
//! \brief Representation of an orthogonal operation on a real-spherical shell -- the spherical counterpart
//! of \c CartesianShellRep.  A real solid harmonic \f$\chi_{l,m}\f$ is a fixed linear combination of the
//! same-degree Cartesian monomials, so the \f$(2l+1)\times(2l+1)\f$ operation rep on the harmonic shell is
//! obtained from the Cartesian rep by projecting through that c2s map.  The harmonic subspace is
//! rotation-invariant, so the projection is EXACT (not a fit):
//! \f[ D_{sph} = (C C^{\mathsf T})^{-1} C\, D_{cart}\, C^{\mathsf T}, \qquad
//!     C = \text{the } n_{sph}\times n_{cart}\text{ coefficient matrix}, \f]
//! with the same transform convention as \c CartesianShellRep: \f$\chi_m(R^{-1}r)=\sum_n D(n,m)\,\chi_n(r)\f$,
//! and \f$R\mapsto D_{sph}(R)\f$ a faithful representation.  Pure math: the c2s coefficients are supplied by
//! the caller (the molecular basis passes its SolidHarmonics), so qcSymmetry stays self-contained + LAPACK-free.
module;
#include <vector>
export module qchem.Symmetry.Molecule.SphericalRep;
export import qchem.Symmetry.Molecule.CartesianRep;   // IVec3(=Monomial), CartesianShellRep, ShellRep; re-exports qcMath Angular

export namespace qchem::Symmetry::Molecule
{

//! The c2s expansion the spherical rep projects through: a list of harmonics, each a list of
//! (Cartesian monomial, coefficient) \c CartTerm's, in the BASIS'S OWN m-ordering (the convention must match
//! the basis these reps adapt -- see the per-basis extractors).  This is the SAME shared type the
//! spherical-Gaussian evaluator produces (\c qchem::Math::SphericalShell), so the extractor hands the basis's
//! own terms across with no re-derivation.  Per-harmonic scale is irrelevant (it cancels in the projection).
using qchem::Math::CartTerm;
using HarmonicC2S = qchem::Math::HarmonicC2S;

//! \brief The operation rep of a real-spherical shell whose harmonics have the Cartesian expansion \a c2s.
//! Built from the Cartesian rep by projecting through the c2s map (the harmonic subspace is
//! rotation-invariant, so the projection is EXACT):
//! \f$\text{Rep}(R) = (C C^{\mathsf T})^{-1} C\, D_{cart}\, C^{\mathsf T}\f$, with \f$C\f$ the
//! \f$n_{sph}\times n_{cart}\f$ coefficient matrix.
class SphericalShellRep : public ShellRep
{
public:
    explicit SphericalShellRep(HarmonicC2S c2s) : itsC2S(std::move(c2s)) {}
    virtual size_t nComponents() const {return itsC2S.size();}
    virtual rmat_t Rep(const rmat3d_t& R) const;
private:
    HarmonicC2S itsC2S;
};

} //namespace

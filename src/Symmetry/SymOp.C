// File: Symmetry/SymOp.C  The structure-NEUTRAL symmetry operation {W|tau} (+ its spin action).
//
// Promoted out of Symmetry/Lattice_3D/Fold.C (R2.17.3, user 2026-08-07).  The struct was always
// generic -- its own doc comment said so -- but its ADDRESS was not: living in the Lattice_3D folder
// made every signature that took one crystal-specific, even where the operation itself is not.  The
// trigger was UnitCell's ops-taking CreateIntegrationMesh: a SITE-ADAPTED quadrature mesh is equally
// meaningful for a molecule, and only the argument's namespace said otherwise.
//
// Placement note (the qcSymmetry taxonomy): this file sits at the ROOT, above the Atom/ Molecule/
// Lattice_3D/ folders, precisely because it belongs to none of them.  That is the same reason Irrep.C
// and Spin.C are here -- the root holds what all three structure families share.
module;
export module qchem.Symmetry.SymOp;
export import qchem.Types;       // rvec3_t
export import qchem.Matrix3D;    // Matrix3D

export namespace qchem::Symmetry
{

//! \brief Spin action carried by a magnetic (Shubnikov) group operation: a plain spatial op
//! leaves the spin channels alone; an op paired with time reversal swaps them (and, on the
//! G-side, conjugates: \f$\tilde\rho(-G)=\tilde\rho(G)^*\f$).  Declared next to the op type, so every
//! consumer is drawn against the spin-native shape from the start (doc/SymmetryUpgradePlan.md §4
//! tier 4a) -- unpol is the \f$\zeta=0\f$ collapse, not the design case.  Purely spatial folds ignore
//! it: \f$\sigma\f$ acts on FIELD channels, never on the point geometry.
//!
//! NON-COLLINEAR heads-up (plan §4/§9): this enum is itself the COLLINEAR collapse of the
//! general spin action.  A spin-space-group op carries a spin ROTATION -- \f$R_s\in SO(3)\f$
//! on the magnetization 3-vector \f$m(r)\f$, equivalently an SU(2) action on the \f$2\times2\f$
//! spinor density -- with time reversal = \f$m\to-m\f$ + conjugation; {None, Flip} is that
//! action restricted to a fixed quantization axis (spirals/canted orderings need the full
//! rotation).  When non-collinear lands, \f$\sigma\f$ grows into that rotation type; consumers
//! must not bake "two scalar channels" into interfaces beyond the collinear tier.
enum class SpinAction { None, Flip };

//! \brief A symmetry operation \f$\{W|\tau\}\f$ -- STRUCTURE-NEUTRAL.
//!
//! \a W is the linear part expressed in the acting space's own basis, and \a tau the translation in
//! the same basis.  For a CRYSTAL that basis is fractional (\a W integer-valued: the direct-space
//! \f$W\f$ for {r} folds, the reciprocal \f$U\f$ for {G}/{k} folds) and \f$\tau\neq0\f$ for a
//! non-symmorphic op -- a screw axis or glide plane, whose translation decides which sites share an
//! orbit, so a bare matrix cannot fold a non-symmorphic raster.  For a MOLECULE or ATOM the basis is
//! Cartesian and \f$\tau=0\f$: a point group fixes a point, so it has no translation part.  A consumer
//! that Cartesianises via \f$A\,W\,A^{-1}\f$ therefore needs no special case -- the molecular \f$A\f$
//! is the identity and the conversion is a no-op.
//!
//! \a tau ACTS on direct-space points (\f$r \to W r + \tau\f$); for {G}/{k} folds it enters only as the
//! edge PHASE \f$e^{+2\pi i\,m\cdot\tau}\f$.
//!
//! \note This is a DISCRETE operation, and a finite list of them is a finite group.  A spherical atom's
//! actual symmetry is the CONTINUOUS group \f$O(3)\f$, which no such list can represent -- and nothing
//! in this code asks it to (user, 2026-08-07).  Where an atom needs a finite op set (a maximally
//! stretched / symmetry-broken configuration, a site group inside a crystal) these ops express it
//! exactly; the continuous case is handled by the \c Symmetry::Atom spherical machinery instead, which
//! works in \f$l,m\f$ rather than by enumerating operations.
struct SymOp
{
    Matrix3D<double> W;                        //!< Linear part, in the acting space's own basis.
    rvec3_t          tau;                      //!< Translation, same basis (zero for a point group).
    SpinAction       sigma = SpinAction::None; //!< Shubnikov spin action (\f$\sigma\f$), §4 tier 4a.
};

} //namespace

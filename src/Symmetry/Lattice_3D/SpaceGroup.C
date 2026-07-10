// File: Symmetry/Lattice_3D/SpaceGroup.C  Crystal space-group detection from a unit cell.
//
// Tier A of the space-group plan (doc/SpaceGroupPlan.md): detect the symmetry operations
// {R|tau} of a crystal from its cell matrix + atom basis, and expose the linear parts as
// the crystal POINT group acting on reciprocal space (k).  Only the point-group part is
// needed for Brillouin-zone reduction, so tau is carried but the k-side consumers ignore it.
//
// Deliberately decoupled from qcStructure (mirroring qchem.Symmetry.Molecule.PointGroup):
// it operates on a cell matrix + a list of fractional AtomSites, so a thin
// UnitCell -> {A, basis} adapter can sit on top without qcSymmetry depending on qcStructure.
module;
#include <vector>
export module qchem.Symmetry.Lattice_3D.SpaceGroup;
export import qchem.Types;       // rvec3_t, ivec3_t, vec3_t
export import qchem.Matrix3D;    // Matrix3D, Vector3D, Determinant, Invert, operators

export namespace qchem::Symmetry::Lattice_3D
{

//! One atom of the crystal basis: a species label (atomic number Z) and its position in
//! FRACTIONAL cell coordinates \f$f\f$ (\f$r = A f\f$).  Two sites are symmetry-interchangeable
//! only if they share a species.
struct AtomSite
{
    int     species;
    rvec3_t f;        //!< Fractional cell coordinates.
};

//! One crystal symmetry operation \f${W|\tau}\f$: an atom at fractional \f$f\f$ maps to
//! \f$W f + \tau\ (\mathrm{mod}\ 1)\f$.  \a W is the linear part expressed in the lattice
//! (fractional-direct) basis -- an integer matrix with \f$\det = \pm1\f$ -- stored as
//! double for convenience.  \a tau is the fractional translation, reduced into \f$[0,1)^3\f$
//! (zero for a symmorphic operation).
struct SpaceGroupOp
{
    Matrix3D<double> W;      //!< Linear part in lattice coordinates (integer-valued).
    rvec3_t          tau;    //!< Fractional translation in \f$[0,1)^3\f$.
};

//! \brief The detected space group of a crystal: its \f${W|\tau}\f$ operations, plus the
//! derived crystal point group acting on \f$k\f$.
//!
//! A space group \f$\approx\f$ point group \f$\ltimes\f$ lattice translations.  For
//! Brillouin-zone reduction (Tier A) only the linear parts \a W matter: in the reciprocal
//! fractional basis \a W acts on \f$k\f$ as \f$(W^{-1})^\top\f$ (so that \f$k'\cdot r' =
//! k\cdot r\f$).  Time reversal \f$k\to-k\f$ is an additional reciprocal-space operation for
//! non-magnetic systems.
class SpaceGroup
{
public:
    //! \brief Detect the space group from a cell matrix and a fractional atom basis.
    //! \param A     Cell matrix; columns are the lattice vectors \f$a_i\f$ (a.u.).
    //! \param basis Atoms in fractional cell coordinates.
    //! \param tol   Absolute tolerance (fractional/length) for matching positions and metrics.
    //!
    //! Enumerates the holohedry (lattice point group: integer ops preserving the metric
    //! \f$M = A^\top A\f$), then keeps those linear parts \a W that admit a fractional
    //! translation \f$\tau\f$ making \f${W|\tau}\f$ a symmetry of the basis.  Assumes a
    //! primitive cell (one \f$\tau\f$ coset per \a W).
    static SpaceGroup Detect(const Matrix3D<double>& A,
                             const std::vector<AtomSite>& basis, double tol=1e-4);

    size_t Order() const {return itsOps.size();}                   //!< Number of \f${W|\tau}\f$ operations.
    const std::vector<SpaceGroupOp>& Ops() const {return itsOps;}  //!< The space-group operations.
    bool   isSymmorphic() const;                                   //!< True if every \f$\tau = 0\f$ (mod lattice).

    //! The crystal POINT group: the linear parts \a W as integer matrices in lattice
    //! (fractional-direct) coordinates.  Order = Order() for a symmorphic group.
    std::vector<Matrix3D<double>> PointGroupOps() const;

    //! \brief The reciprocal-space integer operations acting on \f$k\f$ (fractional reciprocal
    //! coordinates): \f$U = (W^{-1})^\top\f$ for each point-group \a W.
    //! \param includeTimeReversal  Also fold in \f$k\to-k\f$ (adds \f$-U\f$ for each \a U;
    //!        a no-op when the group is centrosymmetric).
    std::vector<Matrix3D<double>> ReciprocalPointOps(bool includeTimeReversal=true) const;

    const Matrix3D<double>& CellMatrix() const {return itsA;}

private:
    SpaceGroup(const Matrix3D<double>& A, std::vector<SpaceGroupOp> ops)
        : itsA(A), itsOps(std::move(ops)) {}

    Matrix3D<double>          itsA;    //!< Cell matrix (columns = lattice vectors).
    std::vector<SpaceGroupOp> itsOps;  //!< The \f${W|\tau}\f$ operations.
};

} // namespace

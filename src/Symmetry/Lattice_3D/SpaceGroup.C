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
export import qchem.Symmetry.Lattice_3D.Fold;   // Fold, SymOp (the shared orbit-fold primitive)

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

//! \brief A reciprocal-space symmetry operation carrying its fractional translation, for the
//! NON-SYMMORPHIC density star-average in G-space (SCATTER form).  \a U is the integer G-index scatter
//! matrix \f$W^\top\f$ (permutes the G-indices exactly, \f$|Um|=|m|\f$); \a tau is the op's fractional
//! translation.  The star-average scatters each input coefficient with the glide phase on its INPUT index
//! (doc/GPWPlan1.md item 5) \f[ \tilde\rho_\mathrm{sym}[U m]\mathrel{+}=\tfrac1{|ops|}\,e^{+2\pi i\,m\cdot\tau}\,\tilde\rho(m), \f]
//! which reproduces the exact projector \f$\tilde\rho_\mathrm{sym}(G)=\tfrac1{|ops|}\sum_{op}e^{+2\pi i(W^{-\top}G)\cdot\tau}\tilde\rho(W^{-\top}G)\f$
//! and reduces to the plain permutation average when every \f$\tau=0\f$ (symmorphic).
struct ReciprocalOp
{
    Matrix3D<double> U;      //!< G-index scatter matrix \f$W^\top\f$ (exact G-index permutation).
    rvec3_t          tau;    //!< Fractional translation of the operation.
};

//! \brief A direct-space symmetry operation carrying its fractional translation, for the
//! NON-SYMMORPHIC density star-average on the real-space raster.  \a W is the integer linear part
//! (real space acts by \f$r\to W r + \tau\f$); the raster gathers from voxel \f$W\cdot g + N\tau\f$
//! (\f$N\tau\f$ = the sublattice shift in voxels), the real-space partner of the \f$e^{+2\pi i(Um)\cdot\tau}\f$
//! G-space phase.
struct DirectOp
{
    Matrix3D<double> W;      //!< Integer linear part in lattice coordinates.
    rvec3_t          tau;    //!< Fractional translation of the operation.
};

//! \brief RUN-LEVEL symmetry-imposition policy (doc/SymmetryUpgradePlan.md §3).  Every symmetry
//! reduction -- the k-fold, the {G}/{r} folds, the stream reduction, the density star-average --
//! IMPOSES the group it folds under: exact iff the density actually HAS that group, silently wrong
//! (too-high-symmetry, too-high-energy) when the ground state wants to break it.  So imposition is
//! an OPT-IN speed contract the caller asserts, never a silent default: out of the box \c None
//! (a free run imposes nothing and the order-parameter diagnostic reports the symmetry the SCF
//! actually found).  ONE policy value feeds every reduction surface of a run; subgroup selection
//! (the SSB-search / ordering-enumeration workflows) extends the enum, not the plumbing.
//!
//! What the caller asserts to (§3): on an imposed run the density defect \f$\|\rho-P_G\rho\|\f$ is
//! \f$\approx 0\f$ BY CONSTRUCTION (every iterate is projected), so the diagnostic cannot audit it
//! -- the RELEASE-CHECK (re-converge under a released group with a symmetry-broken seed) is the
//! audit that certifies no spontaneous symmetry breaking was suppressed.
struct SymmetryPolicy
{
    enum class Impose { None, FullGroup };   //!< Subgroup arrives with the SSB/ordering workflows.
    Impose impose = Impose::None;            //!< Default: impose NOTHING (impose-on-assert).
    bool Imposes() const {return impose != Impose::None;}
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
    //! \param symmorphicOnly  Keep only the linear parts of the \f$\tau=0\f$ (symmorphic) operations.
    //!        This is the SAFE guard for a W-only BZ reduction (density symmetrized without the
    //!        \f$e^{-iG\cdot\tau}\f$ phase): for a symmorphic crystal it changes nothing (every \f$\tau=0\f$),
    //!        while for a non-symmorphic one it drops the glide/screw ops -- folding is merely reduced, never
    //!        wrong (doc/GPWPlan1.md item 3).
    std::vector<Matrix3D<double>> ReciprocalPointOps(bool includeTimeReversal=true,
                                                     bool symmorphicOnly=false) const;

    //! \brief The DIRECT-space integer operations \f$W\f$ (real space acts by \f$r\to W r\f$) -- the reciprocal
    //! partner of \c ReciprocalPointOps, returned straight from \c SpaceGroupOp::W (no invert/transpose).
    //! \param includeTimeReversal  Also add \f$-W\f$ (the real-space partner of \f$k\to-k\f$).  OFF by DEFAULT:
    //!        time reversal is a RECIPROCAL/Bloch concept; in real space \f$-W\f$ imposes a FALSE INVERSION on a
    //!        (real) density, so a caller symmetrizing a real-space raster MUST leave it off -- only the crystal
    //!        point group symmetrizes a real density (doc/GPWPlan1.md item 3).  (In G-space time reversal is
    //!        instead fine, since \f$\tilde\rho(-G)=\tilde\rho(G)^*\f$ -- hence it defaults ON there.)
    //! \param symmorphicOnly  the same \f$\tau=0\f$ W-only guard as \c ReciprocalPointOps.
    std::vector<Matrix3D<double>> DirectPointOps(bool includeTimeReversal=false,
                                                 bool symmorphicOnly=false) const;

    //! \brief The FULL crystal point group as \f${U|\tau}\f$ reciprocal ops (NO symmorphic guard, NO time
    //! reversal) for the non-symmorphic G-space density star-average.  Every \f${W|\tau}\f$ contributes its
    //! reciprocal linear part paired with \f$\tau\f$; the glide phase \f$e^{+2\pi i(Um)\cdot\tau}\f$ makes the
    //! average exact for diamond-type (Fd-3m) crystals, and it collapses to \c ReciprocalPointOps (no TR needed)
    //! when the group is symmorphic.  Time reversal is deliberately absent: the density carries the crystal's
    //! own point symmetry, not an imposed inversion (doc/GPWPlan1.md item 5).
    std::vector<ReciprocalOp> ReciprocalOps() const;

    //! \brief The FULL crystal point group as \f${W|\tau}\f$ direct ops for the non-symmorphic real-space raster
    //! star-average.  The direct partner of \c ReciprocalOps (same op set, linear part \a W kept directly).
    std::vector<DirectOp> DirectOps() const;

    // ------ Orbit folds (doc/SymmetryUpgradePlan.md §2a) -- the ENCAPSULATED fold surface:
    // the caller passes NO ops, the group forwards its own op set to the shared Fold
    // primitive.  (The §3 imposed-subgroup POLICY, when it lands, plugs in here -- one
    // place decides which ops every reduction folds under.)

    //! \brief Fold a Monkhorst-Pack k-grid under the reciprocal ops WITH time reversal (the
    //! k-star).  Raw index = KMesh linear order (ix outer, iz inner).
    Fold FoldKMesh(const ivec3_t& N, const rvec3_t& shift) const;

    //! \brief Fold an explicit G-index list under the full \f$\{U|\tau\}\f$ reciprocal set
    //! (exact integer path, NO time reversal -- the density star-average op set).  Edge op
    //! \a i maps rep \f$m \to U_i m\f$; the member's glide phase for a scatter is
    //! \f$e^{+2\pi i\,m\cdot\tau_i}\f$ (the \c SymmetrizeGMap convention).
    Fold FoldGVectors(const std::vector<ivec3_t>& m) const;

    //! \brief Fold direct-space FRACTIONAL points under the full \f$\{W|\tau\}\f$ set
    //! (tolerance path, NO time reversal -- a real density carries the crystal's own point
    //! symmetry only).  \f$\tau\f$ acts: \f$r \to W r + \tau\f$ (non-symmorphic sites fold).
    Fold FoldPoints(const std::vector<rvec3_t>& pts, double tol = 1e-6) const;

    const Matrix3D<double>& CellMatrix() const {return itsA;}

private:
    SpaceGroup(const Matrix3D<double>& A, std::vector<SpaceGroupOp> ops)
        : itsA(A), itsOps(std::move(ops)) {}

    Matrix3D<double>          itsA;    //!< Cell matrix (columns = lattice vectors).
    std::vector<SpaceGroupOp> itsOps;  //!< The \f${W|\tau}\f$ operations.
};

} // namespace

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
//!
//! \a spin (OPTIONAL -- default 0) is the COLLINEAR magnetic decoration for the Shubnikov
//! classification (\c ShubnikovOps, doc/SymmetryUpgradePlan.md §7 step 7 S1): \f$+1/-1\f$ = the
//! site's moment along/against the shared quantization axis, \f$0\f$ = a non-magnetic site (a
//! closed-shell species, where "flipped" is meaningless and must not constrain the ops).  WHICH
//! sites are magnetic is seed/library knowledge, not geometry -- the caller decorates; plain
//! space-group detection ignores the field entirely.
struct AtomSite
{
    int     species;
    rvec3_t f;        //!< Fractional cell coordinates.
    int     spin = 0; //!< Collinear label: +1/-1 magnetic sublattice, 0 non-magnetic.
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

//! \brief The direct-space face of a reciprocal op: \f$U = W^\top \Rightarrow W = U^\top\f$, same
//! \f$\tau\f$.  Lives HERE so the reciprocal\f$\leftrightarrow\f$direct convention stays with the
//! group -- a consumer holding only the imposed \f$\{U|\tau\}\f$ set recovers the \f$\{W|\tau\}\f$
//! action on real-space points without re-deriving the transpose relation.
inline DirectOp DirectOf(const ReciprocalOp& op) {return {Transpose(op.U), op.tau};}

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

    //! \brief The SITE GROUP (stabilizer) of fractional position \a f: the \f$\{W|\tau\}\f$ ops
    //! with \f$W f + \tau \equiv f \pmod 1\f$.  Always a subgroup; order divides \c Order(), and
    //! \f$|\mathrm{orbit}(f)| = \mathrm{Order}()/|\mathrm{stabilizer}|\f$.  This is the group that
    //! fixes an atom's Becke ANGULAR grid (doc/SymmetryUpgradePlan.md §6 T2: the site group decides
    //! how many and which directions); a generic position returns just the identity.
    std::vector<SpaceGroupOp> SiteStabilizer(const rvec3_t& f, double tol = 1e-6) const;

    //! \brief The LITTLE GROUP of fractional crystal momentum \a kFrac: the \f${W|\tau}\f$ ops whose
    //! reciprocal action fixes \f$k\f$ modulo the reciprocal lattice, \f$(W^{-1})^\top k \equiv k \pmod 1\f$.
    //! Always a subgroup; the full group at \f$\Gamma\f$.  This is the op set a \f$k\f$-resolved reduction
    //! (the T3 stream fold at \f$k\ne\Gamma\f$, doc/SymmetryUpgradePlan.md §6b/T3.4) may fold under: ops
    //! moving \f$k\f$ relate DIFFERENT k-blocks (the IBZ star bookkeeping), not terms within one block.
    std::vector<SpaceGroupOp> LittleGroup(const rvec3_t& kFrac, double tol = 1e-9) const;
    //! \c LittleGroup as \f${W|\tau}\f$ \c DirectOp s -- the currency \c SetStreamSymmetryOps takes.
    std::vector<DirectOp>     LittleGroupDirectOps(const rvec3_t& kFrac, double tol = 1e-9) const;

    //! \brief The COLLINEAR magnetic (Shubnikov) group of a spin-DECORATED basis
    //! (doc/SymmetryUpgradePlan.md §7 step 7, increment S1).  Classify against \a decorated
    //! (same cell; \c AtomSite::spin carries the collinear pattern): an op that maps every site
    //! onto one of EQUAL spin is kept with \f$\sigma=\f$\c None; one that maps every site onto one
    //! of OPPOSITE spin is kept with \f$\sigma=\f$\c Flip (non-magnetic \f$s=0\f$ sites satisfy
    //! both, so an undecorated basis reproduces the grey group, all \c None); an op that does
    //! neither is DROPPED (it is not a symmetry of the magnetic crystal).  \f$M = H + \sigma\cdot
    //! (G\setminus H)\f$ -- always a group; closure is gated in the unit tests.
    //!
    //! COSET COMPLETENESS -- why this does NOT just tag \c Ops(): \c Detect keeps ONE \f$\tau\f$
    //! coset per \f$W\f$ (primitive-cell assumption), but a magnetically DOUBLED cell (MnO AFM-II
    //! = the chemical cell doubled along [111]) has TWO chemical-lattice cosets per \f$W\f$, and
    //! the extra one matters most: the pure ANTI-TRANSLATION \f$\{E|\tfrac12\tfrac12\tfrac12\}
    //! \cdot\f$\c Flip IS the sublattice mirror \f$m_1=-m_2\f$ the imposition exists to enforce.
    //! So this method re-enumerates EVERY valid \f$\tau\f$ per detected \f$W\f$ against the
    //! decorated basis itself; \c Detect and every grey-group consumer stay untouched.
    //!
    //! The returned ops are DIRECT-frame \f$\{W|\tau,\sigma\}\f$ \c SymOp s (the \f$\sigma\f$ slot
    //! reserved by §4 tier 4a).  \f$\sigma\f$ acts on the FIELD channels -- \c Flip swaps
    //! \f$\rho_\uparrow\leftrightarrow\rho_\downarrow\f$ -- never on the point geometry, so the
    //! spatial folds (mesh, {G}, streams) consume \c W|tau unchanged (increment S2/S3 wiring).
    std::vector<SymOp> ShubnikovOps(const std::vector<AtomSite>& decorated, double tol = 1e-4) const;

    const Matrix3D<double>& CellMatrix() const {return itsA;}
    //! The reciprocal cell matrix \f$B = 2\pi A^{-\top}\f$ (columns \f$b_i\f$; \f$G = B\,m\f$).
    Matrix3D<double> ReciprocalMatrix() const;
    //! Cartesian \f$\to\f$ fractional: \f$f = A^{-1} r\f$ (the frame the \f${W|\tau}\f$ ops act in).
    rvec3_t ToFractional(const rvec3_t& r) const {return itsAinv*r;}
    //! Fractional \f$\to\f$ Cartesian: \f$r = A f\f$.
    rvec3_t ToCartesian (const rvec3_t& f) const {return itsA*f;}

private:
    SpaceGroup(const Matrix3D<double>& A, std::vector<SpaceGroupOp> ops)
        : itsA(A), itsAinv(Invert(A)), itsOps(std::move(ops)) {}

    Matrix3D<double>          itsA;    //!< Cell matrix (columns = lattice vectors).
    Matrix3D<double>          itsAinv; //!< \f$A^{-1}\f$ (the fractional-frame converter).
    std::vector<SpaceGroupOp> itsOps;  //!< The \f${W|\tau}\f$ operations.
};

//! \brief The ops COMMON to several candidate magnetic orderings' Shubnikov groups -- the §3
//! imposed-SUBGROUP policy for an ordering SEARCH (doc/SymmetryUpgradePlan.md §3/§7 step 7):
//! impose only what EVERY candidate shares, so the SCF remains free to choose among them.
//! "Shared" means the SAME \f$(W, \tau, \sigma)\f$ TRIPLE: an op that is \c None for one ordering
//! and \c Flip for another (the anti-translation: \c Flip for AFM-II, \c None for FM) is NOT
//! common -- imposing it with either \f$\sigma\f$ would erase one of the candidates.
std::vector<SymOp> CommonOps(const std::vector<std::vector<SymOp>>& groups, double tol = 1e-6);

} // namespace

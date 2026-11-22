// File: Symmetry/Group.C  The Group face: ONE vocabulary for every group the symmetry chain H ⊆ G_ρ ⊆ Go ∩ G_Ven compares.
//
// Design: lane H's Group-face note (side-session notes, 2026-11-19), ruled by the user 2026-11-19 ("1-7) Agree with all of your
// recommendations"; face comments: ONE name `Contains` for a generator and for a group; no public generator list).  This is
// MIGRATION STEP 1: the face, its Generator record, three realizations (finite, crystal, affine stabiliser) built by factories
// from the native kinds, and the exact intersection table.  Nothing else in the tree uses it yet (steps 2-8 move the issuers).
//
// Frame: EVERY operation a Group speaks is CARTESIAN (a.u.): r -> R r + t, R orthogonal.  The native kinds keep their own frames
// (a SpaceGroup's W is fractional) and the factories convert once, R = A W A^-1, t = A tau.
//
// Spaces: a group acts on r, on spin, or jointly on both (magnetic).  A magnetic (Shubnikov) group is NOT a product of an r group
// and a spin group: the SPACES factor (r ⊗ spin) and every element is a pair (g, σ), but the group is a SUBGROUP of G_r × G_spin,
// and the moment pattern selects which pairs survive (note §2.4, user 2026-11-19).  Its elements are ops carrying σ.
module;
#include <memory>
#include <optional>
#include <string>
#include <vector>
export module qchem.Symmetry.Group;
export import qchem.Symmetry.SymOp;               // SymOp {W|tau, sigma}, read here in the CARTESIAN frame
import qchem.Symmetry.Lattice_3D.SpaceGroup;
import qchem.Symmetry.Molecule.PointGroup;
export import qchem.Outcome;

export namespace qchem::Symmetry
{

//! \brief The space(s) a group acts on.  Two groups may form a direct product only across disjoint spaces (user 2026-11-16:
//! "T×O(3) is meaningless (self contradictory) unless they are about different spaces").
enum class Space { R=1, Spin=2, RSpin=3 };

//! \brief Comparison tolerances a group's ISSUER carries (a structure's symmetry tolerance), never a literal at the comparison.
//! Two ops are equal when their linear parts agree to \a angle (max |ΔR_ij|, dimensionless) and their translations to \a length
//! (a.u.).  A comparison between two groups uses the looser of the two (they come from different detections).
struct GroupTolerance
{
    double length;   //!< a.u.: translations, centres, lattice-vector integrality (|A (f - round f)|)
    double angle;    //!< dimensionless: entries of R, and parallelism of axes (|n x n'|)
};

//! \brief The spin group a SpinRotations generator stands for (note §4): SU(2) (restricted, one channel) or U(1) about an axis
//! (collinear, two channels).
struct SpinRotationGroup
{
    enum class Kind { SU2, U1 };
    Kind    kind;
    rvec3_t axis;    //!< the quantization axis (U1 only)
};

//! \brief ONE generator of a group: a discrete operation OR a whole continuous family.  A one-of-N record (ruled 2026-11-19,
//! note Q1): exactly one realization is populated, chosen by the static constructor, and the typed accessors THROW on a
//! mismatch.  Its own clients are the Group realizations; everyone else hands a Generator to \c Group::Contains.
//!
//! A family is exact, never sampled: T(L) is the lattice, O(3)_c and SO(2) about a line are the continuous rotation groups.
class Generator
{
public:
    enum class Kind { Op, Translations, FullRotations, AxialRotations, SpinRotations };

    //! A discrete Cartesian op {R|t, σ}: r -> R r + t, R orthogonal.
    static Generator Op            (const SymOp& g);
    //! T(L): every lattice translation; \a A's COLUMNS are the lattice vectors (a.u.), full rank.
    static Generator Translations  (const Matrix3D<double>& A);
    //! O(3)_c: every orthogonal op (proper AND improper) that fixes the point \a c.
    static Generator FullRotations (const rvec3_t& c);
    //! SO(2) about the line through \a c along \a n: every proper rotation about that axis.  \a n need not be normalized.
    static Generator AxialRotations(const rvec3_t& c, const rvec3_t& n);
    //! The spin group of a spin level (note §4).
    static Generator SpinRotations (const SpinRotationGroup& s);

    Kind GetKind() const {return itsKind;}
    // Typed accessors (the realizations' reads): each THROWS unless its realization is the populated one.
    const SymOp&             GetOp          () const;
    const Matrix3D<double>&  GetLattice     () const;   //!< Translations
    const rvec3_t&           GetCentre      () const;   //!< FullRotations, AxialRotations
    const rvec3_t&           GetAxis        () const;   //!< AxialRotations (unit)
    const SpinRotationGroup& GetSpinGroup   () const;   //!< SpinRotations

private:
    explicit Generator(Kind k) : itsKind(k) {}
    Kind                             itsKind;
    std::optional<SymOp>             itsOp;
    std::optional<Matrix3D<double>>  itsLattice;
    std::optional<rvec3_t>           itsCentre;
    std::optional<rvec3_t>           itsAxis;
    std::optional<SpinRotationGroup> itsSpin;
};

class Group;
using GroupPtr=std::shared_ptr<const Group>;

//! \brief THE FRAME a representation is realized in (Group face step 3b; user 2026-11-21: "We know any radial forces O(3) without creating
//! lists of r,w pairs"): the KIND of a point set or fit fixes the TYPE of group it forces, and the frame supplies only where -- a lattice
//! for a translation family, a centre for a rotation family.  The (r, w) points add nothing.  Either part may be absent: a representation
//! that needs an absent part cannot exist in this frame (a uniform mesh with no lattice, the radial-only set with no centre).
struct Frame
{
    std::optional<Matrix3D<double>> lattice;   //!< columns = the lattice vectors (a.u.)
    std::optional<rvec3_t>          centre;    //!< the one point a rotation family is about
};

//! \brief The FACE every group is asked through (note §2.2).  Data-free; realized by each kind behind the factories below.
//!
//! - \c Contains(Generator) is native and exact: T(L') ⊆ G iff G's translations contain L'; a finite group never contains a
//!   family.
//! - \c Contains(Group) is exact from generators: H ⊆ G iff G contains every generator of H.  No enumeration.
//! - \c Intersect is the one hard request: the side that is finite modulo its translations enumerates, the other filters with
//!   \c Contains; two continuous families meet in a closed table (O(3)_c ∩ O(3)_c' = C∞v about the line c→c'), and two lattices
//!   meet in their common sublattice, found by INTEGER arithmetic (T(L₁) ∩ T(L₂) = T(L₁ ∩ L₂)).  It FAILS (an Outcome with the
//!   reason) for two lattices with no common sublattice: an incommensurate pair is a refusal, never an approximation.
//! - \c ActsOn says which space(s) the group acts on.
//!
//! There is NO public generator list (user 2026-11-19: "What do clients do with std::vector<Generator>?  Why not ask Group to do
//! that operation?").  The one client of a group's generators is \c Contains(Group) itself, below; the intersection table reads
//! each side's canonical form, also protected.  A new client that wants either asks the user first, with its reason.
class Group
{
public:
    virtual ~Group()=default;

    virtual bool  Contains(const Generator& g) const=0;          //!< a discrete op OR a whole family
            bool  Contains(const Group&     h) const;            //!< H ⊆ G (every generator of H is contained)
    Outcome<GroupPtr,std::string> Intersect(const Group& o) const;  //!< G ∩ O; fails only for incommensurate lattices
    virtual Space ActsOn() const=0;
    virtual std::string Describe() const=0;                      //!< one line for a report or a refusal message
    //! \brief Is a LATTICE SUM defined here: does this group hold a full-rank lattice of translations?  What its clients consume:
    //! a periodic seed, the periodic-multipole correctness limit (a lattice sum of site multipoles is conditionally convergent), and a
    //! representation that needs a lattice to exist at all (a uniform cell mesh, the plane-wave ball).  Group face step 3b: replaces
    //! \c Structure::Admits(Translations).
            bool  HasLattice() const;
    //! \brief The frame this group offers a representation (\c Frame): its full translation lattice when it has one (a crystal), and the point
    //! it fixes when it is the rotation group about one point (an atom's O(3)).  A finite group offers neither: its point group fixes no
    //! single site a representation could be pinned at, and it has no lattice.  CLIENT: the pin-30 descriptors (\c qcMesh::Forces,
    //! \c BasisSet::Forces), so a representation that has not been built can say what it would force on this structure.
            Frame MakeFrame() const;

    //! \brief The canonical form the intersection table reads (an implementation detail of this module; see the class note).
    //! One-of-N: a FINITE op list (closed), a CRYSTAL (lattice + one Cartesian coset rep per element mod the lattice), or the
    //! POINTWISE STABILISER of an affine subspace of dimension 0 (O(3)_c) or 1 (C∞v about a line).
    struct Form
    {
        enum class Kind { Finite, Crystal, Stabiliser };
        Kind                 kind;
        std::vector<SymOp>   ops;        //!< Finite: every element; Crystal: one rep per coset of the lattice
        Matrix3D<double>     lattice;    //!< Crystal: columns = the FULL translation lattice
        rvec3_t              centre;     //!< Stabiliser: a point of the fixed subspace
        std::vector<rvec3_t> directions; //!< Stabiliser: 0 (a point) or 1 (a line) unit direction
        GroupTolerance       tol;
    };

protected:
    virtual Form                   GetForm   () const=0;
    //! The generators \c Contains(Group) asks about, derived from the form (Finite: its elements; Crystal: T(L) + the coset reps;
    //! a point: O(3)_c; a line: SO(2) about it + one mirror containing it).
    std::vector<Generator>         Generators() const;
};

// ---------------------------------------------------------------------------------------------------------------------------
//  Factories: the ONLY place the native kinds become Groups (CLAUDE.md: factories concentrate the concrete associations).

//! O(3)_c: the full rotation group about the point \a c (a free atom's spatial group).
GroupPtr PointRotationGroup(const rvec3_t& c, const GroupTolerance& tol);
//! C∞v about the line through \a c along \a n: every orthogonal op fixing the line pointwise.
GroupPtr LineRotationGroup(const rvec3_t& c, const rvec3_t& n, const GroupTolerance& tol);
//! T(L): the pure lattice translations; \a A's columns are the lattice vectors (a.u.).
GroupPtr TranslationGroup(const Matrix3D<double>& A, const GroupTolerance& tol);
//! A finite group from Cartesian ops, CLOSED under composition here (a generating set is enough).  THROWS if the closure does not
//! stay finite (more than \c kMaxFiniteOrder elements: the ops are not a finite group at this tolerance).
GroupPtr FiniteGroup(const std::vector<SymOp>& cartesianOps, const GroupTolerance& tol);
//! A molecular point group from the detected ops about \a origin (each op acts as r -> R (r - origin) + origin).  \a name is what the
//! issuer says the group IS, and \c Describe leads with it: a molecule issuing its abelian subgroup names it so, e.g. "D2h (abelian
//! subgroup of D6h)", so nobody reads it as the full point group (user, 2026-11-20).
GroupPtr FiniteGroup(const std::vector<Molecule::SymOp>& ops, const rvec3_t& origin, const GroupTolerance& tol, const std::string& name={});
//! The space group of a crystal: \a sg's {W|tau} on the cell \a A (columns, a.u.).  A centring translation among the ops (a
//! non-primitive cell) joins the lattice, so the group's translations are its FULL translation subgroup.
GroupPtr CrystalGroup(const Lattice_3D::SpaceGroup& sg, const GroupTolerance& tol);
//! A magnetic (Shubnikov) crystal group: the DIRECT-frame {W|tau, σ} from \c SpaceGroup::ShubnikovOps on the cell \a A.  It acts on
//! (r, spin) jointly when any op carries a flip.
GroupPtr CrystalGroup(const std::vector<SymOp>& fractionalOps, const Matrix3D<double>& A, const GroupTolerance& tol);

//! \brief The tolerance of a group whose generators are EXACT BY CONSTRUCTION, never detected: a fit's or a sampler's \c Forces() is
//! built from its own lattice or centre (Group face step 3), so nothing is measured and roundoff is the only slack.  It matters only
//! when that group is the one asked (\c h.Contains(g)); \c G_struct->Contains(*fit.Forces()) compares at G_struct's own tolerance.
inline constexpr GroupTolerance kExactTolerance{.length=1e-10, .angle=1e-10};
//! \brief The tolerance an ADMISSION question hands a structure's \c MakeSymmetry (Group face step 3b; it was the step-2 bridge's
//! \c theAdmitsBridgeTolerance): pin 30 asks whether G holds T or O(3) about a site, which no realistic geometry decides at the margin, so
//! the structure's usual detection slack (1e-4 bohr, 1e-4 in R) serves every asker.  One constant, so the askers cannot disagree.
inline constexpr GroupTolerance kAdmissionTolerance{.length=1e-4, .angle=1e-4};
//! {e}: the group of a fit or point set that forces nothing beyond its own functions (atom-centred Gaussians, delta functions).
GroupPtr TrivialGroup();

//! The largest finite group \c FiniteGroup will close (O_h × a flip is 96; this leaves room and still stops a runaway closure).
inline constexpr size_t kMaxFiniteOrder=1024;

} //namespace

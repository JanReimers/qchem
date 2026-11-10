// File: Calculation/PotentialEnergySurface.C
//
// qchem::PotentialEnergySurface -- what a geometry optimizer (or an MD / NEB / phonon driver) consumes: the energy, and later
// its gradient, at positions IT chooses, many times, over ONE fixed run setup (doc/CalcNetworkReview.md §3.1, ruling Q2).
// Implemented for molecules (MoleculeSession, qchem.MoleculeCalculation) and for solids (SolidSession, qchem.SolidCalculation): the
// first client that must drive both through one door.
module;
#include <optional>
#include <string>
#include <vector>
export module qchem.PotentialEnergySurface;

export import qchem.Outcome;     // Outcome<T,E>
export import qchem.Types;       // rvec3_t, rmat3d_t
export import qchem.SCFFailure;  // the SCF's own failure, carried by Why::SCF

export namespace qchem
{

//! \brief A GEOMETRY on the surface: what \c Reference() returns and what \c ConvergeAt is given back, in the same units.
struct PESGeometry
{
    //! The CELL: lattice vectors a1, a2, a3 as the COLUMNS of a 3x3 matrix, in BOHR -- the \c UnitCell's own cell matrix
    //! (\c GetCellMatrix).  EMPTY for a molecule.
    std::optional<rmat3d_t> cell;
    std::vector<rvec3_t>    R;   //!< atom positions, CARTESIAN, in BOHR, one per atom in the reference's order
    std::vector<int>        Z;   //!< the atomic number of each atom, in that order (fixed by the setup: a driver reads, never sets, it)
};

//! What a point is asked for besides the energy.  \c stress is the variable-cell seam (not implemented: a request refuses).
struct PESRequest
{
    bool forces = false;
    bool stress = false;
};

//! One converged point: a VALUE, so a driver can keep a history of them.
struct PESPoint
{
    double                              energy = 0.0;   //!< total energy, Ha
    std::optional<std::vector<rvec3_t>> forces;         //!< −dE/dR, Ha/bohr, in the reference's atom order -- present iff requested and supported
    size_t                              iterations = 0; //!< the point's SCF iterations, so a driver can see the carried seed working
};

//! WHY a point produced no energy -- a legitimate outcome the driver must act on (stop, back-track, or ask the user for a new
//! setup), never a crash.  A request that is a PROGRAMMING error (the wrong number of positions) throws instead -- and so does a
//! geometry that breaks an imposed setup's group (lane H, 2026-11-09; OpenWork §2 geometry row): the drivers that move atoms off
//! the group (finite-difference forces, phonons, NEB, the level-2 symmetry probe) run a setup at the lower group, and the level-1
//! optimizer steps through a map that keeps R on the group by construction (\c WyckoffMap, \c CartesianMap given the ops), so an
//! off-group R reaching an imposed session is a wiring defect no driver can act on, not a point to back-track from.
struct PESFailure
{
    enum class Why
    {
        SCF,                 //!< the point's SCF did not converge (details: the SCF's own failure)
        BasisConditioning,   //!< at these positions the frozen basis leaves the conditioning window it was trimmed in (ruling Q3)
        CellChange,          //!< a solid's cell differs from the setup's: the basis and grids cannot follow a cell yet (variable cell)
        NoGradient           //!< forces or stress were asked of a setup whose terms cannot supply them yet
    };
    Why         why = Why::SCF;
    std::string details;
    //! The point's SCF failure itself, when \c why is \c SCF (its \c why, iterations and last-iterate energy, for a driver that
    //! retries with another recipe); empty for the setup refusals.
    std::optional<SCFFailure> scf;
};

//! \brief The energy surface of one fixed run setup, sampled at positions the caller chooses.
//!
//! Data-free and abstract (CLAUDE.md design rule).  An implementer is a SESSION: it owns the current geometry point and the
//! density it carries to the next one, and nothing it does changes the setup -- the basis recipe (its trim frozen), the grids,
//! the pseudopotential set, the imposed group and the +U manifolds are the same at every point.
class PotentialEnergySurface
{
public:
    virtual ~PotentialEnergySurface() = default;
    //! The REFERENCE geometry the setup was resolved at: its cell (a solid's) and positions, units as \c PESGeometry says.
    virtual PESGeometry Reference() const = 0;
    //! \brief Converge at the geometry (\a cell, \a R), seeded from the LAST converged density carried onto it.  Non-const: it
    //! advances the surface's state.
    //! \param cell  lattice vectors as the COLUMNS of a 3x3 matrix, in BOHR; EMPTY for a molecule.  A molecule given a cell, or a
    //!               solid given none, is a broken invariant (throws); a solid's cell that moved is \c PESFailure::CellChange today.
    //! \param R     atom positions, CARTESIAN, in BOHR, \c Reference()'s atom order.  THROWS unless one per atom, and THROWS when an
//!               imposed setup's group does not survive them (a lower group is a NEW setup, ruling Q4).
    virtual Outcome<PESPoint, PESFailure> ConvergeAt(const std::optional<rmat3d_t>& cell, const std::vector<rvec3_t>& R,
                                                     const PESRequest& want = {}) = 0;
};

} // namespace qchem

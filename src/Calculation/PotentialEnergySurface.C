// File: Calculation/PotentialEnergySurface.C
//
// qchem::PotentialEnergySurface -- what a geometry optimizer (or an MD / NEB / phonon driver) consumes: the energy, and later
// its gradient, at positions IT chooses, many times, over ONE fixed run setup (doc/CalcNetworkReview.md §3.1, ruling Q2).
// Implemented for molecules (MoleculeSession, qchem.Calculation) and for solids (SolidSession, qchem.SolidCalculation): the
// first client that must drive both through one door.
module;
#include <optional>
#include <string>
#include <vector>
export module qchem.PotentialEnergySurface;

export import qchem.Structure;   // Structure -- the frame the positions are given in
export import qchem.Outcome;     // Outcome<T,E>
import qchem.Types;              // rvec3_t

export namespace qchem
{

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
//! setup), never a crash.  A request that is a PROGRAMMING error (the wrong number of positions) throws instead.
struct PESFailure
{
    enum class Why
    {
        SCF,                 //!< the point's SCF did not converge (details: the SCF's own failure)
        BasisConditioning,   //!< at these positions the frozen basis leaves the conditioning window it was trimmed in (ruling Q3)
        SymmetryBroken,      //!< the positions break the setup's imposed group (ruling Q4: lowering it is a NEW setup)
        NoGradient           //!< forces or stress were asked of a setup whose terms cannot supply them yet
    };
    Why         why = Why::SCF;
    std::string details;
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
    //! The reference geometry: the species, the atom ORDER the positions follow, and (a solid) the cell, which the setup fixes.
    virtual const Structure& Reference() const = 0;
    //! \brief Converge at \a R (Cartesian bohr, \c Reference()'s atom order), seeded from the LAST converged density carried onto
    //! \a R.  Non-const: it advances the surface's state.  THROWS when \a R does not have one position per atom.
    virtual Outcome<PESPoint, PESFailure> ConvergeAt(const std::vector<rvec3_t>& R, const PESRequest& want = {}) = 0;
};

} // namespace qchem

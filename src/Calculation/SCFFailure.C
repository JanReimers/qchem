// File: Calculation/SCFFailure.C  WHY an SCF attempt has no answer -- the failure half of both facades' Outcome door.
//
// Its own module (2026-10-31, OpenWork "PES / run-setup remainders" (4)): it lived in qchem.SolidCalculation, which the molecular
// facade cannot import, so the molecular Converge stayed a bool.  Both facades -- and the PES face, whose SCF failure carries it --
// now import it from here.
module;
#include <cstddef>
#include <optional>
#include <string>
export module qchem.SCFFailure;
export import qchem.AuditFootprint;   // ResourceRefusal: the numbers a refusal on resources carries

export namespace qchem
{

//! \brief WHY an SCF attempt produced no usable answer (doc/OpenWork.md N1).
//!
//! Carried as a VALUE, not thrown: a non-convergent SCF is a legitimate outcome the caller must act on,
//! not a broken invariant.  \c why exists because the detectors distinguish cases that look identical from
//! the energy alone -- a run that simply ran out of iterations, one whose Hartree term ran away
//! (measured 2026-08-25: Eee 13.5 -> 29.0 Ha), and one that converged to the WRONG STATE, e.g. an imposed
//! AFM cell that relaxed to zero moment, which contradicts its own constraint.
struct SCFFailure
{
    enum class Why
    {
        NotConverged,   //!< ran to the iteration limit with the residual still above tolerance
        ChargeSlosh,    //!< the Hartree term ran away (the low-G charge mode un-damped)
        //! \brief the run did NOT keep the magnetic order it was set up to carry -- either it converged
        //! to zero moment, or the moment rose and then collapsed mid-run.  An imposed Shubnikov group is
        //! an ASSERTION about the answer, so contradicting it is a failed postcondition, not a diagnostic.
        OrderLost,
        //! the FIRST iterate's electronic energy was physically impossible (above the electrons-at-infinity limit 0) and the SCF refused to iterate it: the
        //! Hamiltonian is broken, so no number of iterations would make the answer trustworthy (\c tSCFIterator::UnphysicalStop)
        UnphysicalEnergy,
        //! an audit reference the run needs (the raster pilot's) cannot fit in the memory available: REFUSED ON RESOURCES (CLAUDE.md, three
        //! outcomes).  \c resource carries the numbers and the escape; a refusal mid-pilot keeps the converged run and its record.
        ResourceRefused
    };
    Why         why        = Why::NotConverged;
    size_t      iterations = 0;
    //! The last iterate's energy.  DIAGNOSTIC ONLY -- it is deliberately not reachable as "the energy",
    //! because handing back a plausible number from a failed run is the defect this whole type exists to
    //! prevent (four such numbers were produced in one session: -45.5, -46.3, -56.4, -38.5).
    double      lastEnergy = 0.0;
    std::string details;
    std::optional<ResourceRefusal> resource;   //!< set exactly when \c why is ResourceRefused
};

} // namespace qchem

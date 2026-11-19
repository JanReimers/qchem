// File: Calculation/SetupRefusal.C  WHY a setup refused the options it was given -- the failure half of the setups' Admit doors.
//
// Its own module (2026-11-14, OpenWork "SolidSetup refusals are throws, not Outcomes"), as SCFFailure is: both setups (SolidSetup,
// MoleculeSetup) declare their Admit doors with it, and neither should import the policy layer's helpers (qchem.Policy) to do so.
module;
#include <string>
export module qchem.SetupRefusal;
export import qchem.Outcome;

export namespace qchem
{

//! \brief WHY a setup refused the options it was given (\c SolidSetup::Admit, \c MoleculeSetup::Admit): a VALUE, because a deck, a GUI
//! form or a command line is user input, and a choice the code cannot run is a legitimate failure the caller reports (CLAUDE.md:
//! a call that can legitimately fail returns an Outcome).
struct SetupRefusal
{
    std::string policy;   //!< the policy row's name (the run report's \c policy rows carry it; grep the code for it)
    std::string reason;   //!< the refusal's sentence, as the throwing door would have said it
};

//! \brief WHY a calculation could not be BUILT from an admitted setup (\c AtomCalculation::Prepare, \c MoleculeCalculation::Prepare,
//! \c SolidCalculation::Prepare; V-FACADE-PARITY R9: Setup -> Build -> Go, each phase handing on a proof).  A refusal that fires while the
//! graph is built (an accelerator the cell cannot run, a tensor that will not fit in memory) is a policy's, so it carries the policy row's
//! name beside its sentence, exactly as a \c SetupRefusal does; a broken invariant still throws.
struct BuildRefusal
{
    std::string policy;   //!< the policy row's name (the run report's \c policy rows carry it)
    std::string reason;   //!< the refusal's sentence
};

} //export namespace qchem

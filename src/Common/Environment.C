//! \file Common/Environment.C
//! \brief The RETIRED environment variables (D-ENV step 6): set in a shell, they are ignored and reported.
//!
//! (D-ENV step 4's \c Env(name, legacy) -- a renamed knob read with its old name as a deprecated alias, warning once -- had no caller
//! left once step 6 moved every number into the deck; it was deleted with its warn-once static, D-MUTABLE-STATICS 2026-11-23.)
module;
#include <string>
#include <vector>
export module qchem.Environment;

export namespace qchem
{
//! \brief A RETIRED environment variable (D-ENV step 6): it once set a tier-1/2 value and is now IGNORED -- the deck is the only way to
//! set one.  \c deckKey is the replacement (a dotted deck path, also the \c --set path); EMPTY = removed with no replacement (an A/B switch
//! that is not a run input).
struct RetiredVariable { std::string name; std::string deckKey; };
//! \brief The retired variables that are SET in this process's environment (empty when the environment is clean).  Never honoured.
std::vector<RetiredVariable> RetiredEnvironmentSet();
//! \brief Say ONCE per process, on stderr, which retired variables are set and that they are ignored (and what to use instead).  Called by
//! every run's entry point, so a stale `export` in a shell profile cannot silently change -- or silently fail to change -- a run.
void WarnRetiredEnvironment();
}

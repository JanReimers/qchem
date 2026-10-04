//! \file Common/Environment.C
//! \brief Environment reads with a deprecated alias (D-ENV step 4, 2026-10-03).
//!
//! A setting that is not about the Gaussian-plane-wave basis must not carry the GPW_ prefix (user ruling).  Renaming a
//! knob that banked recipes and run scripts already set must not break them, so every renamed knob is read through
//! \c Env(name, legacy): the new name wins; the old name still works and says ONCE, on stderr, what to use instead.
module;
#include <string>
#include <vector>
export module qchem.Environment;

export namespace qchem
{
//! \brief The value of environment variable \a name, else of its DEPRECATED alias \a legacy (with a one-time stderr notice naming
//! the replacement); nullptr if neither is set.  \a legacy may be null (no alias).
const char* Env(const char* name, const char* legacy=nullptr);

//! \brief A RETIRED environment variable (D-ENV step 6): it once set a tier-1/2 value and is now IGNORED -- the deck is the only way to
//! set one.  \c deckKey is the replacement (a dotted deck path, also the \c --set path).
struct RetiredVariable { std::string name; std::string deckKey; };
//! \brief The retired variables that are SET in this process's environment (empty when the environment is clean).  Never honoured.
std::vector<RetiredVariable> RetiredEnvironmentSet();
//! \brief Say ONCE per process, on stderr, which retired variables are set and that they are ignored (and what to use instead).  Called by
//! every run's entry point, so a stale `export` in a shell profile cannot silently change -- or silently fail to change -- a run.
void WarnRetiredEnvironment();
}

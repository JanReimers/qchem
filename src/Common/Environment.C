// File: Common/Environment.C  Environment reads with a deprecated alias (D-ENV step 4, 2026-10-03).
//
// A setting that is not about the Gaussian-plane-wave basis must not carry the GPW_ prefix (user ruling).  Renaming a
// knob that banked recipes and run scripts already set must not break them, so every renamed knob is read through
// \c Env(name, legacy): the new name wins; the old name still works and says ONCE, on stderr, what to use instead.
module;
export module qchem.Environment;

export namespace qchem
{
//! The value of environment variable \a name, else of its DEPRECATED alias \a legacy (with a one-time stderr notice naming
//! the replacement); nullptr if neither is set.  \a legacy may be null (no alias).
const char* Env(const char* name, const char* legacy=nullptr);
}

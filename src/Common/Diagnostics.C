// File: Common/Diagnostics.C  ONE place for every diagnostic switch (D-ENV step 2, ruled 2026-10-03).
//
// A DIAGNOSTIC prints or records something and NEVER changes a number: a trace, a census, a dump.  Before this module
// each lived as its own `getenv("GPW_...")` in the leaf that used it -- 22 of them, with three different readings of
// "set" (presence, `!=0`, a numeric argument) and no list anywhere, so a mistyped name silently did nothing.  Now:
//
//     QCHEM_DIAGNOSTICS=dm_rank,rss_trace,mesh_ortho=4        one variable, a comma list of `id` or `id=value`
//     QCHEM_DIAGNOSTICS=list                                   print the registry (id, legacy name, what it does)
//
// and each id keeps its OLD environment name as a LEGACY ALIAS (GPW_DM_RANK=1 still works), so no banked recipe breaks.
// An unknown id in QCHEM_DIAGNOSTICS is a WARNING naming the closest ids -- the typo catch the scattered getenvs could
// not give.  A legacy alias is ON when set to anything but "0"; its value (e.g. GPW_MESH_ORTHO=4) is the diagnostic's
// argument.
//
// ⚠ A GLOBAL, ON PURPOSE.  A trace in a leaf function must not need a flag threaded through every constructor, and the
// usual case against globals (hidden coupling in BEHAVIOUR) does not apply: nothing computed may depend on a diagnostic.
// That is the rule that keeps it honest -- if a switch changes a number it is a TYPED OPTION or an A/B hatch, not a
// diagnostic.  Read once at first use; \c Scoped overrides one id for a test (single-threaded; restore on destruction).
module;
#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <vector>
export module qchem.Diagnostics;

export namespace qchem::Diagnostics
{

//! Is diagnostic \a id switched on (QCHEM_DIAGNOSTICS, its legacy environment alias, or a test's \c Scoped)?
//! \a id must be a registered id -- an unregistered one THROWS (a typo in SOURCE is a bug, not a quiet "off").
bool Enabled(const std::string& id);

//! The argument of an enabled diagnostic ("1" when it was just switched on; "4" for \c mesh_ortho=4 or GPW_MESH_ORTHO=4),
//! or nullopt when it is off.
std::optional<std::string> Value(const std::string& id);

//! Every registered diagnostic: id, legacy environment name, one-line description.
struct Entry { std::string id, legacyEnv, what; };
const std::vector<Entry>& Registry();
void Describe(std::ostream& os);

//! The pure parse of a QCHEM_DIAGNOSTICS string: id -> value ("1" when no `=`); ids not in the registry are returned in
//! \a unknown.  A function of its argument so it is unit-testable.
std::map<std::string,std::string> ParseList(const std::string& list, std::vector<std::string>& unknown);

//! TEST ONLY: force \a id on (with \a value) for this object's lifetime, then restore.  Single-threaded.
class Scoped
{
public:
    Scoped(const std::string& id, const std::string& value="1");
    ~Scoped();
    Scoped(const Scoped&)=delete; Scoped& operator=(const Scoped&)=delete;
private:
    std::string itsId;
    std::optional<std::string> itsPrev;
};

} // namespace qchem::Diagnostics

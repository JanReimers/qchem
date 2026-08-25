// File: Common/Outcome.C  The project's vocabulary type for a call that can LEGITIMATELY fail.
//
// THE RULE THIS ENCODES (user, 2026-08-25; doc/OpenWork.md N1).  Two different things get confused under
// the word "error", and they want opposite treatments:
//
//   * A call that can LEGITIMATELY fail on its inputs -- a basis-set lookup that has no entry, a
//     factorisation of a matrix that turns out not to be PSD, an SCF that does not converge.  The caller
//     can and must act on it.  ==> Outcome<T,E>.  No exception; the failure is a VALUE and carries WHY.
//   * A BROKEN INVARIANT -- "this atom-centred mesh has no site blocks", "the composite basis changed
//     mid-run".  No caller can act on it and there is no sensible value to return.  ==> throw.
//
// ⛔ WHY NOT Barton & Nackman's Fallible, which is the obvious prior art.  Its value accessor ASSERTS when
// invalid, and an assert is COMPILED OUT under NDEBUG.  That is exactly the mechanism that kept the
// imposed-mesh site-blocks defect invisible in Release (doc/OpenWork.md item 3) until the assert was
// turned into a throw.  A vocabulary type whose safety evaporates in the build we ship is not safety.
// So: reading the value of a FAILED Outcome THROWS, in every build.
//
// ⚠ AND Outcome IS THE WRONG TOOL FOR A *PHASE*.  "Did the SCF converge?" is ONE fact governing a whole
// FAMILY of queries (energy, terms, density, moments).  Wrapping each accessor in Outcome<T,E> repeats
// that fact N times and asks the caller to check N times.  For a phase, hand back a type that only EXISTS
// on success and put the accessors on IT (SolidCalculation::Converged) -- then there is no check to
// forget, because the accessor is not reachable without the proof.  Outcome carries the WHY of the
// transition; the proof type carries the capability.  (CLAUDE.md: give capabilities only to types that
// have them.)
module;
#include <string>
#include <utility>
#include <variant>
#include <stdexcept>
export module qchem.Outcome;

export namespace qchem
{

//! \brief The result of a call that may legitimately fail: either a \a T or an \a E saying why not.
//!
//! Build with the named factories (never a converting constructor -- \c Outcome<double,double> would be
//! ambiguous, and the ambiguity would be silent).  \c [[nodiscard]] on the type, so DROPPING an Outcome is
//! a compile-time warning at every call site, which is the cheap half of "you must look at this".
template <class T, class E> class [[nodiscard]] Outcome
{
public:
    //! In-place, NEVER via a default-constructed slot: \a T is often a proof type with no default ctor
    //! (that is rather the point of it -- \c SolidCalculation::Converged), and requiring one here would
    //! quietly rule out exactly the types this pairs best with.
    static Outcome Ok  (T v) {return Outcome(std::variant<T,E>(std::in_place_index<0>, std::move(v)));}
    static Outcome Fail(E e) {return Outcome(std::variant<T,E>(std::in_place_index<1>, std::move(e)));}

    bool IsOk() const {return itsSlot.index()==0;}
    explicit operator bool() const {return IsOk();}

    //! The value.  THROWS (not asserts -- see the header) when this Outcome is a failure, so the check
    //! cannot be skipped by building with NDEBUG.
    const T& Value() const
    {
        if (!IsOk()) throw std::runtime_error("Outcome::Value() on a FAILED outcome -- check IsOk() (or "
                                              "operator bool) first, then read Error() for the reason.");
        return std::get<0>(itsSlot);
    }
    T&& TakeValue()
    {
        if (!IsOk()) throw std::runtime_error("Outcome::TakeValue() on a FAILED outcome -- check IsOk() first.");
        return std::move(std::get<0>(itsSlot));
    }
    //! The reason.  Symmetrically, THROWS when this Outcome actually succeeded.
    const E& Error() const
    {
        if (IsOk()) throw std::runtime_error("Outcome::Error() on a SUCCESSFUL outcome -- there is no reason "
                                             "to read; check IsOk() first.");
        return std::get<1>(itsSlot);
    }
    const T& operator* () const {return Value();}
    const T* operator->() const {return &Value();}

private:
    explicit Outcome(std::variant<T,E> slot) : itsSlot(std::move(slot)) {}
    std::variant<T,E> itsSlot;
};

} //export namespace qchem

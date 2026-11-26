// File: Common/tests/Hardening.C  The RelChecked tree's HARDENING is live: an overrun TRAPS there (user ruling 2026-11-26).
//
// RelChecked builds with libc++'s EXTENSIVE hardening and -fsanitize=array-bounds in trap mode (root CMakeLists, QCHEM_RELCHECKED_HARDENING).
// These two death tests plant the overrun classes the flags exist for, one past the end through a RUNTIME index:
//  - a fixed C array INSIDE a class (the 2026-11-03 Hermite3 LMAX overrun: it stayed inside the object, so ASan could not see it);
//  - a std::array member (libc++'s operator[] check).
// Each must die by SIGILL (the trap).  The case is a CANARY: if the flags ever fall off the nightly tree, it fails there.  A tree built
// without the hardening (Release) skips it, because there the overrun silently corrupts the neighbouring member, which is the point.
#include <gtest/gtest.h>
#include <array>
#include <csignal>
#include <version>

namespace
{
#if defined(_LIBCPP_HARDENING_MODE) && defined(_LIBCPP_HARDENING_MODE_EXTENSIVE) && _LIBCPP_HARDENING_MODE == _LIBCPP_HARDENING_MODE_EXTENSIVE
constexpr bool hardened=true;
#else
constexpr bool hardened=false;
#endif

struct CTable   { double c[4]; double after=42; };
struct StdTable { std::array<double,4> c{}; double after=42; };

volatile int onePastTheEnd=4;   // volatile: the index is a runtime value the compiler cannot fold or prove out of range

double WriteCArray()   { CTable t{};  t.c[onePastTheEnd]=7; return t.after; }
double WriteStdArray() { StdTable t;  t.c[onePastTheEnd]=7; return t.after; }
} // namespace

TEST(Hardening, AnOverrunOfACArrayMemberTraps)
{
    if (!hardened) GTEST_SKIP() << "this tree is not hardened (only RelChecked is)";
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    EXPECT_EXIT(WriteCArray(), testing::KilledBySignal(SIGILL), "");
}

TEST(Hardening, AnOverrunOfAStdArrayMemberTraps)
{
    if (!hardened) GTEST_SKIP() << "this tree is not hardened (only RelChecked is)";
    GTEST_FLAG_SET(death_test_style, "threadsafe");
    EXPECT_EXIT(WriteStdArray(), testing::KilledBySignal(SIGILL), "");
}

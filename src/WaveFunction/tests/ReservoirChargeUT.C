// File: ReservoirChargeUT.C  R2.4 (doc/SCFRobustnessPlan.md P2): an occupation reservoir that lost or gained charge throws.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>
import qchem.WaveFunction.Internal.CompositeWF;   // RequireReservoirCharge (unit tests may import internals)

using qchem::WaveFunction::RequireReservoirCharge;

// The measured failure (2026-08-10): the global-mu seed fill of Al left Σ w·n = 2.25 against N = 3, silently in Release.
TEST(ReservoirCharge, R2_4_ALostElectronThrows)
{
    EXPECT_NO_THROW(RequireReservoirCharge(3.0, 3.0, "exact"));
    EXPECT_NO_THROW(RequireReservoirCharge(3.0+1e-12, 3.0, "a Fermi bisection's roundoff"));
    EXPECT_NO_THROW(RequireReservoirCharge(0.0, 0.0, "an empty channel"));
    EXPECT_THROW(RequireReservoirCharge(2.25, 3.0, "the Al seed fill"), std::logic_error);
    EXPECT_THROW(RequireReservoirCharge(7.0, 6.0, "a gained electron"), std::logic_error);
    EXPECT_THROW(RequireReservoirCharge(std::nan(""), 3.0, "a NaN fill"), std::logic_error);
}

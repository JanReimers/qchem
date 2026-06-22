// File: src/Common/tests/Blaze.C  Tests for blazem helpers.
#include "gtest/gtest.h"
import qchem.Blaze;


// VecBuilder: doubling growth + take() shrinks to the logical length, values preserved in order.
TEST(VecBuilder, DoublingAndTake)
{
    blazem::VecBuilder<double> b(2);              // tiny reserve forces several growths
    for (int i = 0; i < 10; ++i) push_back(b, double(i) * 0.5);   // free fn (ADL)
    EXPECT_EQ(b.size(), 10u);

    rvec_t v = b.take();
    ASSERT_EQ(v.size(), 10u);                     // shrunk to logical length, not capacity
    for (int i = 0; i < 10; ++i) EXPECT_DOUBLE_EQ(v[i], double(i) * 0.5);
}

TEST(VecBuilder, Empty)
{
    blazem::VecBuilder<double> b;
    rvec_t v = b.take();
    EXPECT_EQ(v.size(), 0u);
}

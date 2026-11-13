// File: Common/tests/SCFCap.C  The convergence panel's IterationBudget: what it counts, and what it counts APART.
#include "gtest/gtest.h"
import qchem.SCFCap;

using namespace qchem;

TEST(SCFCap, ABudgetCountsTheRunsIterationsAndPrintsAnAuditsApart)
{
    IterationBudget outer(100, "outer (UT)");
    {
        IterationBudget b(10, "the run's own (UT)");
        ReportSCFIterations(7);                     // the run's SCF
        {
            AuditSCF audit;                         // the raster pilot's re-convergence
            ReportSCFIterations(17);
        }
        ReportSCFIterations(2);                     // the run again, after the audit scope closed
        {
            SetupSCF setup;                         // a seed's atomic solve: no budget sees it at all
            ReportSCFIterations(50);
        }
        EXPECT_EQ(b.Used(), 9u)   << "the budget gates the run's own iterations only (TiO2: 43 of 49, then the pilot's 17)";
        EXPECT_EQ(b.Audit(), 17u) << "the audit's iterations are counted apart, and printed";
    }
    EXPECT_EQ(outer.Used(), 9u)   << "an enclosing budget counts what the inner one counts";
    EXPECT_EQ(outer.Audit(), 17u);
}

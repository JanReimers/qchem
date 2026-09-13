// File: src/Hamiltonian/tests/EnergyBreakdown.C  The keyed, role-tagged energy breakdown (V1.12).
#include <gtest/gtest.h>
#include <stdexcept>
import qchem.Energy;

using qchem::EnergyBreakdown;
using qchem::EnergyRole;

// Every contribution sums into the total; the role sums answer the potential / electronic / kinetic questions
// the old fixed members answered; diagnostics are reported and NEVER summed.
TEST(EnergyBreakdown, RoleSumsAndDiagnostics)
{
    EnergyBreakdown e;
    e.Add("Kinetic",  10.0, EnergyRole::Kinetic,   10.0);
    e.Add("Een",     -30.0, EnergyRole::Potential, -30.0);
    e.Add("Eee",       8.0, EnergyRole::Potential,  16.0);   // quadratic: Tr(D J) = 2E
    e.Add("Enn",       3.0, EnergyRole::Constant,    0.0);
    e.Add("MinusTS",  -0.5, EnergyRole::Entropy,     0.0);
    e.AddDiagnostic("EenNL", -7.0);                          // a SUBSET of Een: must not enter any sum

    EXPECT_DOUBLE_EQ(e.GetTotalEnergy(),      10.0-30.0+8.0+3.0-0.5);
    EXPECT_DOUBLE_EQ(e.GetKineticEnergy(),    10.0);
    EXPECT_DOUBLE_EQ(e.GetPotentialEnergy(),  -30.0+8.0+3.0);     // Potential + Constant (the virial's denominator)
    EXPECT_DOUBLE_EQ(e.GetElectronicEnergy(), 10.0-30.0+8.0);     // Kinetic + Potential: no constants, no entropy
    EXPECT_DOUBLE_EQ(e.GetVirial(),           (-30.0+8.0+3.0)/10.0);
    EXPECT_DOUBLE_EQ(e["Een"],   -30.0);
    EXPECT_DOUBLE_EQ(e["Absent"],  0.0);
    EXPECT_DOUBLE_EQ(e.Diagnostic("EenNL"), -7.0);
    EXPECT_FALSE(e.Has("EenNL"));                            // a diagnostic is not a contribution
}

// Add merges by NAME: the two spin channels (or per-irrep pieces) of one term land in one entry, and the
// expectations add too.  op+= is the same merge over a whole breakdown.
TEST(EnergyBreakdown, MergeByName)
{
    EnergyBreakdown up, dn;
    up.Add("Exc", -2.0, EnergyRole::Potential, -2.5);
    dn.Add("Exc", -1.0, EnergyRole::Potential, -1.5);
    up+=dn;
    EXPECT_EQ(up.Terms().size(), 1u);
    EXPECT_DOUBLE_EQ(up["Exc"], -3.0);
    EXPECT_DOUBLE_EQ(*up.Terms()[0].second.TrDV, -4.0);
    EXPECT_THROW(up.Add("Exc", 1.0, EnergyRole::Kinetic), std::logic_error);   // a name has ONE role
}

// The band form: sum f eps + sum (E - Tr(D V)).  The kinetic correction is identically 0 (E == Tr(D T)), so
// for a linear+quadratic breakdown the band total equals the direct one when sum f eps = sum Tr(D V).
TEST(EnergyBreakdown, BandFormMatchesDirectAndThrowsWhenUnclaimed)
{
    EnergyBreakdown e;
    e.Add("Kinetic",  10.0, EnergyRole::Kinetic,   10.0);
    e.Add("Een",     -30.0, EnergyRole::Potential, -30.0);
    e.Add("Eee",       8.0, EnergyRole::Potential,  16.0);
    e.Add("Enn",       3.0, EnergyRole::Constant,    0.0);
    const double sumFEps = 10.0-30.0+16.0;                   // = sum Tr(D V) over the Fock's terms
    EXPECT_DOUBLE_EQ(e.GetBandEnergy(sumFEps), e.GetTotalEnergy());

    e.Add("Exc", -1.0, EnergyRole::Potential);               // no Tr(D V) claimed (a mixed-density XC term)
    EXPECT_THROW(e.GetBandEnergy(sumFEps), std::logic_error);
    EXPECT_DOUBLE_EQ(e.GetTotalEnergy(), 10.0-30.0+8.0+3.0-1.0);   // the direct total is unaffected
}

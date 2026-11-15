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
    up.Add("Eee", -2.0, EnergyRole::Potential, -2.5);
    dn.Add("Eee", -1.0, EnergyRole::Potential, -1.5);
    up+=dn;
    EXPECT_EQ(up.Terms().size(), 1u);
    EXPECT_DOUBLE_EQ(up["Eee"], -3.0);
    EXPECT_DOUBLE_EQ(*up.Terms()[0].second.TrDV, -4.0);
    EXPECT_THROW(up.Add("Eee", 1.0, EnergyRole::Kinetic), std::logic_error);   // a name has ONE role
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

    e.AddXC(-1.0, std::nullopt, "DM", true);                 // no Tr(D V) claimed (a mixed-density XC term)
    EXPECT_THROW(e.GetBandEnergy(sumFEps), std::logic_error);
    EXPECT_DOUBLE_EQ(e.GetTotalEnergy(), 10.0-30.0+8.0+3.0-1.0);   // the direct total is unaffected
}

// THE PIN-28 CERTIFICATE FAILS CLOSED (user 2026-11-15): an XC energy enters only through AddXC, with its certificate; a breakdown with an
// XC energy is certified only if EVERY XC contribution certified rho[D], on either side of a merge; with no XC energy nothing is claimed.
TEST(EnergyBreakdown, TheXCCertificateFailsClosed)
{
    EnergyBreakdown none;
    none.Add("Kinetic", 1.0, EnergyRole::Kinetic, 1.0);
    EXPECT_FALSE(none.charge.xcBuilt);
    EXPECT_THROW(none.Add("Exc", -1.0, EnergyRole::Potential), std::logic_error) << "the XC name has one door, AddXC";

    EnergyBreakdown good, bad;
    good.AddXC(-1.0, std::nullopt, "DM", true);
    EXPECT_TRUE(good.charge.xcBuilt && good.charge.xcRhoCertified);
    bad.AddXC(-0.5, std::nullopt, "Pρ", false);
    EXPECT_TRUE(bad.charge.xcBuilt && !bad.charge.xcRhoCertified);

    EnergyBreakdown a=good; a+=bad;
    EXPECT_FALSE(a.charge.xcRhoCertified) << "one uncertified XC term uncertifies the sum";
    EXPECT_DOUBLE_EQ(a["Exc"], -1.5);
    EnergyBreakdown b=bad;  b+=good;
    EXPECT_FALSE(b.charge.xcRhoCertified) << "in either order";
    EnergyBreakdown c=none; c+=good;
    EXPECT_TRUE(c.charge.xcBuilt && c.charge.xcRhoCertified) << "a side with no XC energy claims nothing";
    good.AddXC(-0.1, std::nullopt, "fitted", false);
    EXPECT_FALSE(good.charge.xcRhoCertified) << "a second XC term in one breakdown must certify too";
}

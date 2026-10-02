// File: src/WaveFunction/tests/EigenRowsUT.C  The level-table ROWS (D10).  Each test below is a bug that was
// previously found only by reading a run log (CleanupHistory2 D10, 2026-08-08/10) -- now a property of a value.
#include "gtest/gtest.h"
#include <memory>

import qchem.WaveFunction.EigenRows;
import qchem.Symmetry.Atom.Internal.SphericalQNs;   // Yl: a constructible atomic symmetry label
import qchem.Types;

using namespace qchem;
using namespace qchem::WaveFunction;
using qchem::Orbitals::EnergyLevel;
using qchem::Orbitals::EnergyLevels;

namespace
{
sym_t Sym(size_t l) { return std::make_shared<const Symmetry::Atom::Internal::Yl>(l); }
EnergyLevel Level(double e, double occ, int degen, size_t n, Spin ms, size_t l) { return EnergyLevel(e, occ, degen, Orbital_QNs(n, ms, Sym(l))); }
EnergyLevels Levels(std::initializer_list<EnergyLevel> ls) { EnergyLevels els; for (const auto& l : ls) els.insert(l); return els; }
}

TEST(EigenRows, OccupationCellKeepsAFractionalOccupationFractional)
{
    EXPECT_EQ(OccupationCell(2.0, 1), "2/1");
    EXPECT_EQ(OccupationCell(0.0, 3), "0/3");
    EXPECT_EQ(OccupationCell(0.996, 1), "1.00/1") << "a smeared 0.996 must show its decimals, not round to an integer configuration";
    EXPECT_EQ(OccupationCell(0.004, 1), "0.00/1") << "...nor 0.004 to \"0/1\"";
    EXPECT_NE(OccupationCell(0.996, 1), "1/1");
}

// The 2026-08-08 bug: the unpolarized table stopped at the first empty level, which under MOM hides a hole BELOW an
// occupied level.  It must run to the highest occupied level and only then stop.
TEST(EigenRows, UnpolarizedTableRunsPastAHoleToTheHighestOccupiedLevel)
{
    const EnergyLevels els=Levels({ Level(-1.29, 0.0, 1, 1, Spin::None, 0),    // an EMPTY level below...
                                    Level(-0.50, 2.0, 1, 2, Spin::None, 0),    // ...an occupied one
                                    Level( 0.75, 0.0, 1, 3, Spin::None, 0),    // a virtual past the frontier
                                    Level( 0.90, 0.0, 1, 4, Spin::None, 0) });
    const auto rows=UnpolarizedEigenRows(els);
    ASSERT_EQ(rows.size(), 2u) << "hole + HOMO kept, the first virtual past the frontier ends the table";
    EXPECT_EQ(rows[0].occ, "0/1");
    EXPECT_EQ(rows[1].occ, "2/1");
}

TEST(EigenRows, UnpolarizedTableStopsByOccupationNotByEnergySign)
{
    // A metal's occupied levels can be POSITIVE (arbitrary energy zero): the old `e>0 => stop` hid all of them.
    const EnergyLevels els=Levels({ Level(0.10, 2.0, 1, 1, Spin::None, 0), Level(0.30, 2.0, 1, 2, Spin::None, 0),
                                    Level(0.50, 0.0, 1, 3, Spin::None, 0) });
    EXPECT_EQ(UnpolarizedEigenRows(els).size(), 2u);
}

// The 2026-08-10 bug (MnO run 29): an ABSENT channel's energy was filled from the other channel, so a row read as a
// level empty at an energy where the opposite spin is occupied, with the splitting printing exactly 0.00000000.
TEST(EigenRows, AnAbsentChannelPrintsDashesNeverTheOtherChannelsNumber)
{
    const EnergyLevels up=Levels({ Level(-0.60, 1.0, 1, 1, Spin::Up, 0) });
    const EnergyLevels dn=Levels({});                                          // the down channel has no such level
    const EnergyLevels combined=Levels({ Level(-0.60, 1.0, 1, 1, Spin::None, 0) });
    const auto rows=PolarizedEigenRows(combined, up, dn);
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0].eUp, "-0.60000000");
    EXPECT_EQ(rows[0].occDn, "--");
    EXPECT_EQ(rows[0].eDn,   "--") << "not the up energy";
    EXPECT_EQ(rows[0].dE,    "--") << "no splitting is computed from a level that does not exist";
    EXPECT_TRUE(rows[0].dnDim);
}

TEST(EigenRows, ExchangeSplittingIsUpMinusDown)
{
    const EnergyLevels up=Levels({ Level(-0.60, 1.0, 1, 1, Spin::Up,   0) });
    const EnergyLevels dn=Levels({ Level(-0.45, 0.0, 1, 1, Spin::Down, 0) });
    const EnergyLevels combined=Levels({ Level(-0.60, 1.0, 1, 1, Spin::None, 0), Level(-0.45, 0.0, 1, 1, Spin::None, 0) });
    const auto rows=PolarizedEigenRows(combined, up, dn);
    ASSERT_EQ(rows.size(), 1u) << "the two combined entries are ONE (n,sym) row";
    EXPECT_EQ(rows[0].dE, "-0.15000000");
    EXPECT_EQ(rows[0].occUp, "1/1");
    EXPECT_EQ(rows[0].occDn, "0/1");
    EXPECT_TRUE(rows[0].dnDim);
}

// A doubly-empty level BELOW the frontier is a HOLE and must be kept; one ABOVE it is a virtual and is dropped.
TEST(EigenRows, ADoublyEmptyLevelIsKeptBelowTheFrontierAndDroppedAbove)
{
    const EnergyLevels up=Levels({ Level(-0.80, 0.0, 1, 1, Spin::Up, 0), Level(-0.60, 1.0, 1, 2, Spin::Up, 0), Level(0.40, 0.0, 1, 3, Spin::Up, 0) });
    const EnergyLevels dn=Levels({ Level(-0.80, 0.0, 1, 1, Spin::Down, 0), Level(-0.55, 1.0, 1, 2, Spin::Down, 0), Level(0.40, 0.0, 1, 3, Spin::Down, 0) });
    const EnergyLevels combined=Levels({ Level(-0.80, 0.0, 2, 1, Spin::None, 0), Level(-0.60, 1.0, 1, 2, Spin::None, 0),
                                         Level(-0.55, 1.0, 1, 2, Spin::None, 0), Level(0.40, 0.0, 2, 3, Spin::None, 0) });
    const auto rows=PolarizedEigenRows(combined, up, dn);
    ASSERT_EQ(rows.size(), 2u) << "the hole at -0.80 and the occupied level; the virtual at +0.40 is dropped";
    EXPECT_EQ(rows[0].occUp, "0/1");
    EXPECT_EQ(rows[1].occUp, "1/1");
}

TEST(EigenRows, SmearedOccupationsStayFractionalInThePolarizedTableToo)
{
    const EnergyLevels up=Levels({ Level(-0.5, 0.996, 1, 1, Spin::Up,   0) });
    const EnergyLevels dn=Levels({ Level(-0.4, 0.004, 1, 1, Spin::Down, 0) });
    const EnergyLevels combined=Levels({ Level(-0.5, 0.996, 1, 1, Spin::None, 0), Level(-0.4, 0.004, 1, 1, Spin::None, 0) });
    const auto rows=PolarizedEigenRows(combined, up, dn);
    ASSERT_EQ(rows.size(), 1u);
    EXPECT_EQ(rows[0].occUp, "1.00/1");
    EXPECT_EQ(rows[0].occDn, "0.00/1");
}

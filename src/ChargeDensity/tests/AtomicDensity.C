// File: UnitTests/AtomicDensityUT.C  Tests for the SAD atomic-density database reader + interpolation.
#include "gtest/gtest.h"
#include <unistd.h>    // getpid: per-PROCESS temp names (ctest runs each case as its own process, in parallel)
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <cctype>
#include <algorithm>
#include <cmath>
#include "nlohmann/json.hpp"   // the seed gates read the library and NIST's tables directly

import qchem.ChargeDensity.AtomicDensity;   // GetAtomicDensity, RadialDensity, RecentredAtomicDensity
import qchem.Types;                          // rvec3_t
using namespace qchem;

using namespace qchem::ChargeDensity;

// The stored grid integrates 4*pi*int r^2 rho dr to the element's electron count (the generator's charge
// check), and the reader's own trapezoid reproduces it.
TEST(AtomicDensity, ChargeIntegratesToNelec)
{
    for (int Z : {1,2,6,8,10})
    {
        RadialDensity rad = GetAtomicDensity(Z, "LDA");
        // Relative bound: the log-radial quadrature has a small (~0.04%) error from the sharp core, which
        // grows with Z in absolute terms -- a Z-scaled tolerance is the honest "integrates to ~Nelec" check.
        EXPECT_NEAR(rad.Charge(), double(Z), Z*2e-3) << "Z=" << Z;
    }
}

// rho(r) is positive, peaks at the core, and decays outward; node interpolation stays in range.
TEST(AtomicDensity, RadialShape)
{
    RadialDensity rad = GetAtomicDensity(8, "LDA");   // oxygen
    EXPECT_GT(rad(0.0), 0.0);
    EXPECT_GT(rad(0.1), rad(1.0));
    EXPECT_GT(rad(1.0), rad(5.0));
    EXPECT_NEAR(rad(1e6), 0.0, 1e-30);                // clamps to 0 far outside the grid
}

// Recentring: rho(|r-R|) -- value at the nucleus equals rad(0); an offset of d equals rad(d).
TEST(AtomicDensity, RecentredMatchesRadial)
{
    auto rad = std::make_shared<const RadialDensity>(GetAtomicDensity(8, "LDA"));
    rvec3_t R(0.3, -1.2, 2.0);
    RecentredAtomicDensity cd(rad, R);
    EXPECT_DOUBLE_EQ(cd(R), (*rad)(0.0));
    const double d = 0.75;
    EXPECT_DOUBLE_EQ(cd(R + rvec3_t(0,0,d)), (*rad)(d));
    EXPECT_DOUBLE_EQ(cd(R + rvec3_t(d,0,0)), (*rad)(d));
}

TEST(AtomicDensity, MissingElementThrows)
{
    EXPECT_THROW(GetAtomicDensity(123, "LDA"), std::runtime_error);
    EXPECT_THROW(GetAtomicDensity(8, "NoSuchFunctional"), std::runtime_error);
}

// ---------------- spin-resolved pairs (doc/SCFSeedingPlan.md sec 10, the spin-SAD tables) ----------------

static constexpr const char* vdb = "atomic_valence_densities.json";

// The Hund entries carry an up-majority pair whose sum is the stored rho and whose difference integrates
// to the Hund moment 2S: Mn q7 neutral (4s^2 3d^5, 2S=5), Mn2+ (3d^5, 2S=5 -- the physical TM-cation
// removal order, 4s first), O neutral (2p^4 triplet, 2S=2).
TEST(AtomicDensity, SpinPairChargesAndMoments)
{
    struct Case { int Z, Nelec; double moment; };
    for (Case c : {Case{25,7,5.0}, Case{25,5,5.0}, Case{8,6,2.0}})
    {
        ASSERT_TRUE(HasAtomicSpinPair(c.Z, "LDA", vdb, c.Nelec)) << "Z="<<c.Z<<" Nelec="<<c.Nelec;
        auto [up,dn] = GetAtomicSpinPair(c.Z, "LDA", vdb, c.Nelec);
        // Nelec-scaled bound, same reason as ChargeIntegratesToNelec: the reader's log-mesh quadrature
        // carries a small sharp-core error on top of the generator's own ~1e-3.
        EXPECT_NEAR(up.Charge()+dn.Charge(), double(c.Nelec), c.Nelec*1e-3) << "Z="<<c.Z; // sum = charge state
        EXPECT_NEAR(up.Charge()-dn.Charge(), c.moment,        c.Nelec*1e-3) << "Z="<<c.Z; // difference = 2S
        EXPECT_GE(up.Charge(), dn.Charge());                                        // up-majority convention

        // The pair sums to the spin-summed rho the unpolarized readers see (same entry, same grid).  The
        // three arrays are emitted independently at 12 significant digits, so consistency is RELATIVE.
        RadialDensity tot = GetAtomicDensity(c.Z, "LDA", vdb, c.Nelec);
        for (double r : {0.05, 0.3, 1.0, 3.0})
            EXPECT_NEAR(up(r)+dn(r), tot(r), 1e-9*tot(r)+1e-15) << "Z="<<c.Z<<" r="<<r;
    }
}

// Closed-shell (or spin-agnostic legacy) entries carry no pair: Has is false and Get throws.
TEST(AtomicDensity, SpinPairAbsentForClosedShells)
{
    EXPECT_FALSE(HasAtomicSpinPair( 8, "LDA", vdb, 8));   // O2- (closed p^6, generated unpolarized)
    EXPECT_FALSE(HasAtomicSpinPair(14, "LDA", vdb, 4));   // Si (legacy spin-agnostic entry)
    EXPECT_THROW(GetAtomicSpinPair(14, "LDA", vdb, 4), std::runtime_error);
}

// D-SEED1: the seed library is keyed by (Z, functional) but a PSEUDOPOTENTIAL VARIANT (q) is part of a seed density's
// identity (Li q1 and q3 are different densities).  A fixture library (an absolute path replaces the data dir) with
// Li q1's neutral (Nelec 1), Li q3's neutral (Nelec 3) and Li q3's +2 ion (Nelec 1): "the neutral Li" is ambiguous
// and so is "the Nelec=1 Li" -- both must THROW rather than pick the first in file order; "Nelec=3" is unique.
TEST(AtomicDensity, AmbiguousPseudopotentialVariantThrowsInsteadOfPickingTheFirst)
{
    const std::string path=(std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_seed_ambiguity.json")).string();
    {
        auto entry=[](int nelec, int q)
        {
            return "{\"Z\":3,\"functional\":\"LDA\",\"Nelec\":"+std::to_string(nelec)+",\"q\":"+std::to_string(q)
                 + ",\"grid\":{\"rmin\":1e-3,\"rmax\":10.0},\"rho\":[0.1,0.05,0.01]}";
        };
        std::ofstream f(path);
        f<<"["<<entry(1,1)<<","<<entry(3,3)<<","<<entry(1,3)<<"]";
    }
    EXPECT_THROW(GetAtomicDensity(3, "LDA", path),    std::runtime_error) << "two neutrals (q1, q3)";
    EXPECT_THROW(GetAtomicDensity(3, "LDA", path, 1), std::runtime_error) << "q1's neutral vs q3's +2 ion, both Nelec=1";
    EXPECT_NO_THROW(GetAtomicDensity(3, "LDA", path, 3)) << "Nelec=3 is unique";
    EXPECT_FALSE(HasAtomicDensity(3, "LDA", path, 4)) << "no such state is absence, not an error";
}

// C4 step 3: the same fixture, now asked WITH the PP variant q (SeedKey::Q) -- the ambiguities dissolve and q1 / q3
// pick different entries (different rho), and a q the library lacks is absence.
TEST(AtomicDensity, TheQKeyPicksTheVariantsOwnEntry)
{
    const std::string path=(std::filesystem::path(testing::TempDir())/(std::to_string(::getpid())+"_seed_qkey.json")).string();
    {
        auto entry=[](int nelec, int q, double rho0)
        {
            return "{\"Z\":3,\"functional\":\"LDA\",\"Nelec\":"+std::to_string(nelec)+",\"q\":"+std::to_string(q)
                 + ",\"grid\":{\"rmin\":1e-3,\"rmax\":10.0},\"rho\":["+std::to_string(rho0)+",0.05,0.01]}";
        };
        std::ofstream f(path);
        f<<"["<<entry(1,1,0.1)<<","<<entry(3,3,0.3)<<","<<entry(1,3,0.7)<<"]";
    }
    const double r1=GetAtomicDensity(3,"LDA",path,-1,1)(1e-3), r3=GetAtomicDensity(3,"LDA",path,-1,3)(1e-3);
    EXPECT_NEAR(r1,0.1,1e-12) << "q1's neutral";
    EXPECT_NEAR(r3,0.3,1e-12) << "q3's neutral (Nelec==q), not the q3 +2 ion";
    EXPECT_NEAR(GetAtomicDensity(3,"LDA",path,1,3)(1e-3),0.7,1e-12) << "q3's +2 ion (Nelec 1) is addressable once q says which Nelec=1";
    EXPECT_FALSE(HasAtomicDensity(3,"LDA",path,-1,5)) << "a q the library lacks is absence";
    const SeedKey key("LDA",{{3,1}});
    EXPECT_EQ(key.Q(3),1);
    EXPECT_EQ(key.Q(8),-1);
}

// ---------------- THE ALL-ELECTRON LIBRARY, every entry judged (OpenWork §3 "The default SAD seed stops at Ne", lane B 2026-11-13) ----------
//
// scripts/gen-ae-seeds writes atomic_densities.json from the atom facade and stamps each generated entry with its Etot, its configuration
// (NIST's, dftdata/configurations) and the oracle it was judged on.  These gates re-judge the FILE, so a hand edit or a regeneration
// that drifts fails here: the charge of every entry, the spin pair of every polarized one, and the energy of every stamped one against
// NIST's own tables, read independently of the generator.

namespace
{
nlohmann::json ReadJson(const std::filesystem::path& p)
{
    std::ifstream f(p);
    if (!f) throw std::runtime_error("cannot read "+p.string());
    nlohmann::json j; f>>j; return j;
}
const nlohmann::json& AELibrary()
{
    static const nlohmann::json lib=ReadJson(std::filesystem::path(CHARGEDENSITY_DATA_PATH)/"atomic_densities.json");
    return lib;
}
//! Hund's 2S of a configuration text ("[Ar] 3d^5 4s^1"): each partly filled subshell of 2l+1 orbitals holds min(n, 2(2l+1)-n) unpaired
//! electrons (the noble-gas core is closed, so only the explicit subshells can be open).
int HundUnpaired(const std::string& config)
{
    int unpaired=0;
    for (size_t i=0;i<config.size();i++)
    {
        const char l=config[i];
        if (!(l=='s'||l=='p'||l=='d'||l=='f') || i==0 || !std::isdigit((unsigned char)config[i-1]) || i+2>=config.size() || config[i+1]!='^') continue;
        const int deg=2*(2*(std::string("spdf").find(l))+1);
        const int n=std::stoi(config.substr(i+2));
        unpaired+=std::min(n, deg-n);
    }
    return unpaired;
}
} // namespace

// Coverage: every element NIST tabulates (Z 1-92) has a neutral all-electron LDA seed, and the reader returns it at its own charge.
TEST(AtomicDensity, EveryElementThroughUraniumHasANeutralAllElectronSeed)
{
    for (int Z=1; Z<=92; Z++)
    {
        ASSERT_TRUE(HasAtomicDensity(Z, "LDA")) << "Z=" << Z;
        EXPECT_NEAR(GetAtomicDensity(Z, "LDA").Charge(), double(Z), Z*2e-3) << "Z=" << Z;
    }
}

// Every entry: the stored charge is its Nelec (the log grid's quadrature error is ~4e-4 relative); every POLARIZED entry: the pair is
// up-majority, sums to rho, and its moment is the Hund 2S of the configuration it states.
TEST(AtomicDensity, EveryAllElectronEntryHoldsItsChargeAndItsSpinPairItsHundMoment)
{
    for (const auto& e : AELibrary())
    {
        const int Z=e.at("Z"), Nelec=e.at("Nelec");
        EXPECT_NEAR(e.at("charge").get<double>(), double(Nelec), Nelec*1e-3) << "Z=" << Z;
        if (!e.contains("rho_up")) continue;
        const auto& rho=e.at("rho"); const auto& up=e.at("rho_up"); const auto& dn=e.at("rho_dn");
        ASSERT_EQ(up.size(), rho.size()) << "Z=" << Z;
        ASSERT_EQ(dn.size(), rho.size()) << "Z=" << Z;
        for (size_t i=0;i<rho.size();i++)
            EXPECT_NEAR(up[i].get<double>()+dn[i].get<double>(), rho[i].get<double>(), 1e-9*rho[i].get<double>()+1e-15) << "Z=" << Z << " node " << i;
        ASSERT_TRUE(e.contains("config")) << "Z=" << Z << ": a polarized entry states the configuration its moment is judged on";
        const double m=e.at("moment");
        EXPECT_GE(m, 0.0) << "Z=" << Z << ": up-majority";
        EXPECT_NEAR(m, double(HundUnpaired(e.at("config"))), 1e-2) << "Z=" << Z << " " << e.at("config").get<std::string>();
        auto [ru,rd]=GetAtomicSpinPair(Z, "LDA", "atomic_densities.json", Nelec);
        EXPECT_NEAR(ru.Charge()-rd.Charge(), m, Nelec*1e-3) << "Z=" << Z << ": the reader sees the same pair";
    }
}

// Every STAMPED entry (Etot): its energy is NIST's for the same functional and spin treatment -- LSD for a polarized entry, LDA for an
// unpolarized one -- to 1e-5 relative.  This is what catches a generator that converged to another configuration, which the charge cannot.
TEST(AtomicDensity, EveryStampedAllElectronEntryIsNISTsAtomToTenPpm)
{
    const auto common=std::filesystem::path(COMMON_DATA_PATH);
    const nlohmann::json lda=ReadJson(common/"nistLDA.json"), lsd=ReadJson(common/"nistLSD.json");
    for (const auto& e : AELibrary())
    {
        if (!e.contains("Etot")) continue;
        const int Z=e.at("Z");
        const bool pol=e.contains("rho_up");
        const double E=e.at("Etot"), Eref=(pol ? lsd : lda).at(std::to_string(Z)).at("Etot");
        EXPECT_LE(std::abs(E-Eref), 1e-5*std::abs(Eref)) << "Z=" << Z << (pol ? " vs NIST LSD " : " vs NIST LDA ") << Eref << ", E " << E;
        EXPECT_TRUE(e.value("converged", false)) << "Z=" << Z;
    }
}

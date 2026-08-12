// File: src/ChargeDensity/tests/JointPulay.C  ONE Pulay history over BOTH spin channels.
//
// A memoryless FILTER (Kerker's G^2/(G^2+G0^2)) is linear and channel-diagonal, so running one per spin
// channel is identical to running it jointly.  An EXTRAPOLATOR is not: its coefficients come from a
// least-squares fit over a residual HISTORY, so two independent fits give c_up != c_dn and the extrapolated
// pair (Sum c_up_i rho_up_i, Sum c_dn_i rho_dn_i) is a state that NEVER OCCURRED on the trajectory -- each
// channel still conserves its own charge, but the MOMENT comes out an arbitrary synthesised combination of
// the history moments.  That is the MnO AFM-II ejection (doc/SymmetryUpgradePlan.md sec 7 step 7): a run
// sitting at E ~ -61.25, m_stag ~ 0.35 thrown to -49 with the moment collapsed, then climbing back.
//
// The cure is one B summed over the channels and ONE coefficient vector -- spin is just another irrep,
// exactly as the Fock-space DIIS already sums its B over spatial irreps.  These tests drive PulayMixer
// through both paths on a hand-built trajectory and check it against an INDEPENDENT reference written here
// (its own map arithmetic, its own DIIS calls), so a change to the mixer's algebra cannot silently drag the
// oracle with it.
//
// The mixers are built with alpha=1 and G0=0, i.e. the Kerker filter is identically 1 and the step is a pure
// passthrough: the mixed state IS rho~_out*, so every assertion reads the EXTRAPOLATION directly with no
// filter algebra in the way.
#include "gtest/gtest.h"
#include <complex>
#include <memory>
#include <vector>

import qchem.ChargeDensity.DensityMixer;   // PulayMixer, KerkerMixer, GField, MixJointly
import qchem.ChargeDensity.FourierMixCD;   // FourierMixCD::RhoTilde
import qchem.UnitCell;                     // UnitCell + MakeReciprocalCell
import qchem.ReciprocalLattice;            // ReciprocalLattice
import qchem.Math.DIIS;                    // the bordered solve (the reference's own calls)
import qchem.Blaze;                        // rsmat_t / rvec_t / blazem::zero
import qchem.Types;                        // dcmplx, ivec3_t

using namespace qchem;
using namespace qchem::ChargeDensity;

namespace
{

ReciprocalLattice Recip(double a) { UnitCell cell(a); return ReciprocalLattice(cell.MakeReciprocalCell()); }

const ivec3_t Z(0,0,0), A(1,0,0), B(0,1,0), C(1,1,0);

ΔG_Map Field(double z, double a, double b, double c)
{
    ΔG_Map m;
    m[Z]=dcmplx(z,0.0); m[A]=dcmplx(a,0.0); m[B]=dcmplx(b,0.0); m[C]=dcmplx(c,0.0);
    return m;
}

// --- the reference's own map arithmetic (deliberately NOT the module's MapSub/MapCombine/MapInnerRe) ---
ΔG_Map Sub(const ΔG_Map& x, const ΔG_Map& y)
{
    ΔG_Map r=x;
    for (const auto& [k,v]:y) r[k]-=v;
    return r;
}
ΔG_Map Combine(const std::vector<ΔG_Map>& maps, const rvec_t& c)
{
    ΔG_Map r;
    for (size_t i=0;i<maps.size();++i) for (const auto& [k,v]:maps[i]) r[k]+=c[i]*v;
    return r;
}
double Inner(const ΔG_Map& x, const ΔG_Map& y)      // Re Sum_{G!=0} conj(x)*y -- G=0 carries the charge
{
    double s=0.0;
    for (const auto& [k,v]:x)
        if (!(k==Z)) if (auto it=y.find(k); it!=y.end()) s+=std::real(std::conj(v)*it->second);
    return s;
}
rsmat_t Overlap(const std::vector<ΔG_Map>& res)
{
    const size_t n=res.size();
    rsmat_t Bm=blazem::zero<double>(n);
    for (size_t i=0;i<n;++i) for (size_t j=i;j<n;++j) Bm(i,j)=Inner(res[i],res[j]);
    return Bm;
}
double MaxDiff(const ΔG_Map& x, const ΔG_Map& y)
{
    double d=0.0;
    for (const auto& [k,v]:x) d=std::max(d, std::abs(v - (y.count(k) ? y.at(k) : dcmplx(0.0))));
    for (const auto& [k,v]:y) if (!x.count(k)) d=std::max(d,std::abs(v));
    return d;
}

//! THE TRAJECTORY: the freshly collocated rho~_out of each channel, per iteration.  Deliberately asymmetric
//! (the up channel's residual lives mostly on G=A, the down channel's on G=B) -- that asymmetry is exactly
//! what makes two independent fits disagree.
struct Trajectory
{
    std::vector<ΔG_Map> up{Field(0.520,0.20,0.14,0.030), Field(0.505,0.26,0.11,0.045), Field(0.500,0.22,0.13,0.035)};
    std::vector<ΔG_Map> dn{Field(0.510,0.16,0.22,0.040), Field(0.505,0.12,0.27,0.025), Field(0.500,0.14,0.25,0.030)};
    ΔG_Map seedUp=Field(0.50,0.30,0.10,0.05);
    ΔG_Map seedDn=Field(0.50,0.10,0.30,0.02);
    size_t Steps() const { return up.size(); }
};

//! The reference mixer: alpha=1, G0=0 (so mixed == rho~_out*), history depth 8, no priming.  \a joint
//! selects the ONE-B-summed-over-channels solve; otherwise each channel solves its own B (the defect).
struct Reference
{
    std::vector<ΔG_Map> mixedUp, mixedDn;   //!< the mixed state after each step
    std::vector<rvec_t> cUp, cDn;           //!< the coefficients used at each step (empty while n<2)

    Reference(const Trajectory& t, bool joint)
    {
        std::vector<ΔG_Map> outsU, outsD, resU, resD;
        ΔG_Map inU=t.seedUp, inD=t.seedDn;
        for (size_t k=0;k<t.Steps();++k)
        {
            outsU.push_back(t.up[k]); resU.push_back(Sub(t.up[k],inU));
            outsD.push_back(t.dn[k]); resD.push_back(Sub(t.dn[k],inD));
            rvec_t cu, cd;
            if (resU.size()>=2)
            {
                if (joint)
                {
                    rsmat_t Bj=Overlap(resU); Bj+=Overlap(resD);      // B_ij = Sum_sigma <r_i,r_j>_sigma
                    cu=cd=Math::DIIS::Coefficients(Math::DIIS::Bordered(Bj));
                }
                else
                {
                    cu=Math::DIIS::Coefficients(Math::DIIS::Bordered(Overlap(resU)));
                    cd=Math::DIIS::Coefficients(Math::DIIS::Bordered(Overlap(resD)));
                }
                inU=Combine(outsU,cu); inD=Combine(outsD,cd);          // alpha=1 => mixed == rho~_out*
            }
            else { inU=t.up[k]; inD=t.dn[k]; }                         // n<2 => the plain un-extrapolated step
            mixedUp.push_back(inU); mixedDn.push_back(inD);
            cUp.push_back(cu); cDn.push_back(cd);
        }
    }
};

std::unique_ptr<PulayMixer> MakePulay(const ReciprocalLattice& r, const ΔG_Map& seed, double charge)
{   // alpha=1, G0=0 (filter == 1), depth 8, start 0 (no priming); no fit basis => no raw-raster shadow.
    return std::make_unique<PulayMixer>(1.0, 0.0, 8, 0, nullptr, r, seed, charge, rvec_t{});
}

} //namespace

// The CURE: staged through MixJointly, the two channels are extrapolated with ONE coefficient vector solved
// from the SUMMED B -- so the mixed pair is Sum c_i (rho_up_i, rho_dn_i), a point in the affine span of the
// PAIRS the trajectory actually visited, and its moment is Sum c_i m_i with the same weights.
TEST(JointPulay, OneHistoryOverBothChannels)
{
    const ReciprocalLattice recip=Recip(10.0);
    const Trajectory t;
    const Reference  ref(t,/*joint*/true);
    auto up=MakePulay(recip,t.seedUp,8.0), dn=MakePulay(recip,t.seedDn,6.0);

    for (size_t k=0;k<t.Steps();++k)
    {
        MixJointly({up.get(),dn.get()}, {GField{t.up[k],rvec_t{}}, GField{t.dn[k],rvec_t{}}});
        EXPECT_LT(MaxDiff(up->Mixed().RhoTilde(), ref.mixedUp[k]), 1e-12) << "up channel, step " << k;
        EXPECT_LT(MaxDiff(dn->Mixed().RhoTilde(), ref.mixedDn[k]), 1e-12) << "dn channel, step " << k;
    }
    // The moment of the extrapolated state is the SAME combination of history moments as the density.
    const rvec_t& c=ref.cUp.back();
    ASSERT_EQ(c.size(), t.Steps());
    std::vector<ΔG_Map> moments;
    for (size_t i=0;i<t.Steps();++i) moments.push_back(Sub(t.up[i],t.dn[i]));
    EXPECT_LT(MaxDiff(Sub(up->Mixed().RhoTilde(), dn->Mixed().RhoTilde()), Combine(moments,c)), 1e-12);
}

// THE DEFECT, kept measurable: one PulayMixer per channel driven independently (MixField, the single-channel
// path) fits two DIFFERENT coefficient vectors, and the moment it synthesises is off the trajectory.  The
// first half of this test is also the regression that the single-channel path itself is unchanged.
TEST(JointPulay, SplitHistoryFitsDifferentCoefficientsAndMovesTheMoment)
{
    const ReciprocalLattice recip=Recip(10.0);
    const Trajectory t;
    const Reference  split(t,/*joint*/false), joint(t,/*joint*/true);
    auto up=MakePulay(recip,t.seedUp,8.0), dn=MakePulay(recip,t.seedDn,6.0);

    for (size_t k=0;k<t.Steps();++k)
    {
        up->MixField(GField{t.up[k],rvec_t{}});          // the single-channel path: stage, solve MY B, apply
        dn->MixField(GField{t.dn[k],rvec_t{}});
        EXPECT_LT(MaxDiff(up->Mixed().RhoTilde(), split.mixedUp[k]), 1e-12) << "up channel, step " << k;
        EXPECT_LT(MaxDiff(dn->Mixed().RhoTilde(), split.mixedDn[k]), 1e-12) << "dn channel, step " << k;
    }
    // The two fits genuinely disagree on this trajectory...
    const rvec_t &cu=split.cUp.back(), &cd=split.cDn.back();
    ASSERT_EQ(cu.size(), t.Steps());
    double dc=0.0; for (size_t i=0;i<cu.size();++i) dc=std::max(dc,std::fabs(cu[i]-cd[i]));
    EXPECT_GT(dc, 1e-3) << "the trajectory is not asymmetric enough to exercise the defect";
    // ...and the moment they synthesise is NOT the one a single fit gives -- the ejection, in miniature.
    const ΔG_Map mSplit=Sub(up->Mixed().RhoTilde(), dn->Mixed().RhoTilde());
    const ΔG_Map mJoint=Sub(joint.mixedUp.back(), joint.mixedDn.back());
    EXPECT_GT(MaxDiff(mSplit,mJoint), 1e-4);
}

// The property the composite branches on: an extrapolator ANSWERS with its staging face, a memoryless filter
// does not.  This is what makes "one leaf per channel" legitimate for Kerker and illegitimate for Pulay.
TEST(JointPulay, OnlyExtrapolatorsCarryHistory)
{
    const ReciprocalLattice recip=Recip(10.0);
    const ΔG_Map seed=Field(0.5,0.3,0.1,0.05);
    auto pulay=MakePulay(recip,seed,8.0);
    KerkerMixer kerker(0.4, 1.0, nullptr, std::make_shared<FourierMixCD>(seed,recip,8.0), rvec_t{});
    EXPECT_NE(pulay->History(),  nullptr);
    EXPECT_EQ(kerker.History(),  nullptr);
}

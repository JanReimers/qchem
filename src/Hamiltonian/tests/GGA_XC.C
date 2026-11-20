// File: Hamiltonian/tests/GGA_XC.C  The GGA functionals from libxc (GGA stage 1, user ruling 2026-11-20): libxc's derivatives against a
// finite difference of libxc's own energy (trust but verify), and the batch face against the point face for the LDAs.
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>
#include "gtest/gtest.h"

import qchem.Hamiltonian.Internal.ExFunctional;      // XCInput / XCOutput, the batch face
import qchem.Hamiltonian.Internal.Libxc_GGA;         // the GGA under test
import qchem.Hamiltonian.Internal.Libxc_LDA;
import qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.VWN_Correlation;
import qchem.Symmetry.Spin;
import qchem.Types;
using namespace qchem;
using namespace qchem::Hamiltonian;

namespace
{
//! One point of (ρ↑, ρ↓, σ↑↑, σ↑↓, σ↓↓).
struct P { double up, dn, uu, ud, dd; };

//! The test points: densities over five decades (the low-density tail included), spin polarizations from ζ=0 to ζ=0.8 of either
//! sign, and reduced gradients s = |∇ρ_σ|/(2ρ_σ)^{4/3} from 0.05 to 3 (well past PBE's enhancement knee), with the two channel gradients
//! at 60° (σ↑↓ = |∇ρ↑||∇ρ↓| cos 60°).
std::vector<P> Points()
{
    std::vector<P> ps;
    for (double rho : {1e-4, 3e-3, 0.05, 0.7, 12.0})
        for (double zeta : {0.0, 0.3, -0.55, 0.8})
            for (double s : {0.05, 0.4, 1.2, 3.0})
            {
                const double up=0.5*rho*(1+zeta), dn=0.5*rho*(1-zeta);
                const double gu=s*std::pow(2*up,4.0/3.0), gd=s*std::pow(2*dn,4.0/3.0);
                ps.push_back({up, dn, gu*gu, 0.5*gu*gd, gd*gd});
            }
    return ps;
}

XCInput Input(const std::vector<P>& ps)
{
    const size_t n=ps.size();
    XCInput in{rvec_t(n), rvec_t(n), rvec_t(n), rvec_t(n), rvec_t(n)};
    for (size_t g=0;g<n;++g) { in.up[g]=ps[g].up; in.dn[g]=ps[g].dn; in.sigUU[g]=ps[g].uu; in.sigUD[g]=ps[g].ud; in.sigDD[g]=ps[g].dd; }
    return in;
}

XCOutput Eval(const ExFunctional& f, const XCInput& in)
{
    XCOutput out; out.Reset(in); f.Evaluate(in, out); return out;
}

//! e at one point with one argument displaced by \a h.
double EAt(const ExFunctional& f, P p, int arg, double h)
{
    double* a[5]={&p.up, &p.dn, &p.uu, &p.ud, &p.dd};
    *a[arg]+=h;
    return Eval(f, Input({p})).e[0];
}
} // namespace

// THE DERIVATIVES libxc hands us (v_ρσ = ∂e/∂ρ_σ, v_σστ = ∂e/∂σ_στ) are its energy's derivatives: a four-point central difference of e in
// each of the five arguments, step relative to the argument, at every test point, for every GGA part the models use.  MEASURED 2026-11-20
// (libxc d00e5d0e): worst |v - fd|/scale = 2.6e-11 over the set (PBEsol correlation; PBE exchange 3.1e-12), the difference's own O(h⁴)
// truncation plus O(ε/h) roundoff.  The gate is 1e-9: 40x that, tight enough that a factor-of-two slip in any v (the σ↑↓ convention, a
// missing ρ weight) fails by eight orders.
TEST(LibxcGGA, EachDerivativeIsItsEnergysFiniteDifference)
{
    for (int id : {101, 130, 116, 133, 106, 131})   // PBE x/c, PBEsol x/c, B88, LYP
    {
        const Libxc_GGA f(id);
        const std::vector<P> ps=Points();
        const XCOutput out=Eval(f, Input(ps));
        for (size_t g=0;g<ps.size();++g)
        {
            const double v[5]={out.vUp[g], out.vDn[g], out.vsUU[g], out.vsUD[g], out.vsDD[g]};
            const double x[5]={ps[g].up, ps[g].dn, ps[g].uu, ps[g].ud, ps[g].dd};
            for (int a=0;a<5;++a)
            {
                const double h=1e-4*std::abs(x[a]);
                if (!(h>0.0)) continue;
                const double fd=(8.0*(EAt(f,ps[g],a,h)-EAt(f,ps[g],a,-h))-(EAt(f,ps[g],a,2*h)-EAt(f,ps[g],a,-2*h)))/(12.0*h);
                // The scale: this derivative, or (for a near-zero one) e over its argument -- a relative test where the derivative is large,
                // an absolute one on the energy's own scale where it vanishes (exchange's v_σ↑↓ is exactly 0).
                const double scale=std::max(std::abs(v[a]), std::abs(out.e[g]/x[a]));
                EXPECT_NEAR(v[a], fd, 1e-9*scale) << "libxc id " << id << ", point " << g << ", argument " << a;
            }
        }
    }
}

// THE UNPOLARIZED COLLAPSE the term will feed (ρ↑ = ρ↓ = ρ/2, σ_στ = σ/4): both channels answer alike, as they must for a ζ = 0 run.
TEST(LibxcGGA, TheUnpolarizedCollapseAnswersBothChannelsAlike)
{
    for (int id : {101, 130, 106, 131})
    {
        const Libxc_GGA f(id);
        std::vector<P> ps;
        for (double rho : {1e-3, 0.2, 5.0}) for (double sig : {1e-6, 0.03, 2.0}) ps.push_back({rho/2, rho/2, sig/4, sig/4, sig/4});
        const XCOutput out=Eval(f, Input(ps));
        for (size_t g=0;g<ps.size();++g)
        {
            EXPECT_NEAR(out.vUp[g], out.vDn[g], 1e-13*std::abs(out.vUp[g])) << id;
            EXPECT_NEAR(out.vsUU[g], out.vsDD[g], 1e-13*std::abs(out.vsUU[g])+1e-300) << id;
        }
    }
}

// THE BATCH FACE OF AN LDA IS ITS POINT FACE, number for number, PART BY PART: the default loops the spin-native face, so every LDA term
// moved onto the batch later keeps each functional's numbers.  The SUM associates differently (the composite's point face is
// up·Σε↑ + dn·Σε↓, the batch Σ(up·ε↑ + dn·ε↓)): its e agrees to roundoff, its v (a plain sum either way) bit for bit.  And the gradient
// question: no LDA reads ∇ρ; any sum holding a GGA does.
TEST(LibxcGGA, AnLDAsBatchFaceIsItsPointFaceAndOnlyAGGAReadsTheGradient)
{
    const std::vector<std::shared_ptr<ExFunctional>> parts{std::make_shared<SlaterExchange>(2.0/3.0), std::make_shared<VWN_Correlation>(),
                                                           std::make_shared<Libxc_LDA>(20)};
    const CompositeExFunctional lda(parts);
    XCInput in;
    in.up=rvec_t{0.3, 1e-3, 2.0, 0.4}; in.dn=rvec_t{0.3, 4e-3, 0.5, 0.0};   // ζ = 0, ζ < 0, ζ > 0, fully polarized
    for (const auto& f : parts)
    {
        XCOutput out; out.Reset(in); f->Evaluate(in, out);
        EXPECT_EQ(out.vsUU.size(), 0u) << "an LDA input carries no sigma, and the output no v_sigma";
        for (size_t g=0;g<in.size();++g)
        {
            EXPECT_EQ(out.e  [g], f->GetExcDensity(in.up[g], in.dn[g]));
            EXPECT_EQ(out.vUp[g], f->GetVxc(in.up[g], in.dn[g], Spin::Up));
            EXPECT_EQ(out.vDn[g], f->GetVxc(in.up[g], in.dn[g], Spin::Down));
        }
    }
    XCOutput out; out.Reset(in); lda.Evaluate(in, out);
    for (size_t g=0;g<in.size();++g)
    {
        const double e=lda.GetExcDensity(in.up[g], in.dn[g]);
        EXPECT_NEAR(out.e[g], e, 4e-16*std::abs(e)) << "the sum's e: roundoff of the association only";
        EXPECT_EQ(out.vUp[g], lda.GetVxc(in.up[g], in.dn[g], Spin::Up));
        EXPECT_EQ(out.vDn[g], lda.GetVxc(in.up[g], in.dn[g], Spin::Down));
    }
    EXPECT_FALSE(lda.ConsumesDensityGradient());
    const CompositeExFunctional pbe({std::make_shared<Libxc_GGA>(101), std::make_shared<Libxc_GGA>(130)});
    EXPECT_TRUE(pbe.ConsumesDensityGradient());
    const CompositeExFunctional mixed({std::make_shared<SlaterExchange>(2.0/3.0), std::make_shared<Libxc_GGA>(131)});
    EXPECT_TRUE(mixed.ConsumesDensityGradient()) << "one GGA part makes the sum read the gradient";
}

// THE REFUSALS: a GGA has no ρ-only point value (a route that reaches one should have been refused by policy), an input without σ is a
// broken invariant, and an id that is not a pure GGA (an LDA, a hybrid with exact exchange) is not this evaluator's.
TEST(LibxcGGA, APointFaceAMissingGradientAndANonGGAIdThrow)
{
    const Libxc_GGA f(101);
    EXPECT_THROW((void)f.GetVxc(0.1), std::logic_error);
    EXPECT_THROW((void)f.GetEpsXc(0.1, 0.05, Spin::Up), std::logic_error);
    XCInput lda; lda.up=rvec_t{0.1}; lda.dn=rvec_t{0.1};
    XCOutput out; out.Reset(lda);
    EXPECT_THROW(f.Evaluate(lda, out), std::logic_error);
    EXPECT_THROW(Libxc_GGA(1),   std::invalid_argument) << "LDA_X is an LDA";
    EXPECT_THROW(Libxc_GGA(402), std::invalid_argument) << "B3LYP is a hybrid";
    EXPECT_THROW(Libxc_GGA(-7),  std::invalid_argument) << "no such id";
}

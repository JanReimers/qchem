// File: Hamiltonian/tests/GGA_XC.C  The GGA functionals from libxc (GGA stage 1, user ruling 2026-11-20): libxc's derivatives against a
// finite difference of libxc's own energy (trust but verify), the LDA batch face against the LDA point face, and the gradient as a TYPE
// (user 2026-11-20): a GGA has no rho-only face and its input cannot be made without sigma -- checked at COMPILE time.
#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include "gtest/gtest.h"

import qchem.Hamiltonian.Internal.ExFunctional;      // LDAInput / LDAOutput, the LDA batch face
import qchem.Hamiltonian.Internal.GGAFunctional;     // GGAInput / GGAOutput
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

GGAInput Input(const std::vector<P>& ps)
{
    const size_t n=ps.size();
    rvec_t up(n), dn(n), uu(n), ud(n), dd(n);
    for (size_t g=0;g<n;++g) { up[g]=ps[g].up; dn[g]=ps[g].dn; uu[g]=ps[g].uu; ud[g]=ps[g].ud; dd[g]=ps[g].dd; }
    return GGAInput(up, dn, uu, ud, dd);
}

GGAOutput Eval(const GGAFunctional& f, const GGAInput& in)
{
    GGAOutput out(in); f.Evaluate(in, out); return out;
}

//! e at one point with one argument displaced by \a h.
double EAt(const GGAFunctional& f, P p, int arg, double h)
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
        const GGAOutput out=Eval(f, Input(ps));
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
        const GGAOutput out=Eval(f, Input(ps));
        for (size_t g=0;g<ps.size();++g)
        {
            EXPECT_NEAR(out.vUp[g], out.vDn[g], 1e-13*std::abs(out.vUp[g])) << id;
            EXPECT_NEAR(out.vsUU[g], out.vsDD[g], 1e-13*std::abs(out.vsUU[g])+1e-300) << id;
        }
    }
}

// THE LDA BATCH FACE IS ITS POINT FACE, number for number, PART BY PART: the default loops the spin-native face, so every LDA term
// moved onto the batch later keeps each functional's numbers.  The SUM associates differently (the composite's point face is
// up·Σε↑ + dn·Σε↓, the batch Σ(up·ε↑ + dn·ε↓)): its e agrees to roundoff, its v (a plain sum either way) bit for bit.
TEST(LibxcGGA, AnLDAsBatchFaceIsItsPointFace)
{
    const std::vector<std::shared_ptr<ExFunctional>> parts{std::make_shared<SlaterExchange>(2.0/3.0), std::make_shared<VWN_Correlation>(),
                                                           std::make_shared<Libxc_LDA>(20)};
    const CompositeExFunctional lda(parts);
    const LDAInput in{rvec_t{0.3, 1e-3, 2.0, 0.4}, rvec_t{0.3, 4e-3, 0.5, 0.0}};   // ζ = 0, ζ < 0, ζ > 0, fully polarized
    for (const auto& f : parts)
    {
        LDAOutput out(in); f->Evaluate(in, out);
        for (size_t g=0;g<in.size();++g)
        {
            EXPECT_EQ(out.e  [g], f->GetExcDensity(in.up[g], in.dn[g]));
            EXPECT_EQ(out.vUp[g], f->GetVxc(in.up[g], in.dn[g], Spin::Up));
            EXPECT_EQ(out.vDn[g], f->GetVxc(in.up[g], in.dn[g], Spin::Down));
        }
    }
    LDAOutput out(in); lda.Evaluate(in, out);
    for (size_t g=0;g<in.size();++g)
    {
        const double e=lda.GetExcDensity(in.up[g], in.dn[g]);
        EXPECT_NEAR(out.e[g], e, 4e-16*std::abs(e)) << "the sum's e: roundoff of the association only";
        EXPECT_EQ(out.vUp[g], lda.GetVxc(in.up[g], in.dn[g], Spin::Up));
        EXPECT_EQ(out.vDn[g], lda.GetVxc(in.up[g], in.dn[g], Spin::Down));
    }
}

// THE GRADIENT IS A TYPE (user 2026-11-20: "If b the compile fails"): these are COMPILE-TIME facts, so a regression fails the build.
// A GGA is not an ExFunctional, so it has no rho-only (or rho_up, rho_dn) face to call; it takes no LDA input; and a GGA input cannot be
// made from densities alone.  Nothing here is a run-time refusal to test, because the wrong call does not exist.
template <class F>            concept HasRhoOnlyFace   = requires(const F& f) { f.GetVxc(0.1); };
template <class F>            concept HasSpinPointFace = requires(const F& f) { f.GetExcDensity(0.1, 0.1); };
template <class F, class In, class Out> concept TakesInput = requires(const F& f, const In& in, Out& out) { f.Evaluate(in, out); };
static_assert(!std::is_base_of_v<ExFunctional, Libxc_GGA>,             "a GGA is not an LDA-faced functional");
static_assert(!std::is_base_of_v<ExFunctional, CompositeGGAFunctional>);
static_assert( HasRhoOnlyFace<SlaterExchange> && HasSpinPointFace<SlaterExchange>, "(the probes see an LDA's faces)");
static_assert(!HasRhoOnlyFace<Libxc_GGA>,                              "no rho-only face on a GGA");
static_assert(!HasSpinPointFace<Libxc_GGA>,                            "no (rho_up, rho_dn) face on a GGA");
static_assert( TakesInput<Libxc_GGA, GGAInput, GGAOutput> && TakesInput<SlaterExchange, LDAInput, LDAOutput>, "(the probes see the right inputs)");
static_assert(!TakesInput<Libxc_GGA, LDAInput, LDAOutput>,                      "a GGA takes no LDA input");
static_assert(!TakesInput<SlaterExchange, GGAInput, GGAOutput>,                 "an LDA takes no GGA input");
static_assert(!std::is_constructible_v<GGAInput, rvec_t, rvec_t>,      "sigma is REQUIRED: no GGA input from densities alone");
static_assert(!std::is_default_constructible_v<GGAInput>);

// WHAT IS LEFT AT RUN TIME is what only a run-time value can get wrong: a libxc id (user-reachable through XCFunctional) that is not a
// pure GGA, and a GGA input whose five arrays disagree in length.
TEST(LibxcGGA, ANonGGAIdAndAMisSizedInputThrow)
{
    EXPECT_THROW(Libxc_GGA(1),   std::invalid_argument) << "LDA_X is an LDA";
    EXPECT_THROW(Libxc_GGA(402), std::invalid_argument) << "B3LYP is a hybrid";
    EXPECT_THROW(Libxc_GGA(-7),  std::invalid_argument) << "no such id";
    EXPECT_THROW(GGAInput(rvec_t{0.1}, rvec_t{0.1}, rvec_t{0.0}, rvec_t{}, rvec_t{0.0}), std::invalid_argument) << "sigma_ud has no entry";
    EXPECT_THROW(CompositeGGAFunctional({}), std::invalid_argument);
}

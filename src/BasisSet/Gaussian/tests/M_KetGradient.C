// File: BasisSet/Gaussian/tests/M_KetGradient.C  THE 2-CENTRE KERNELS' KET-CENTRE DERIVATIVE, Enzyme's (GaussianRF::KetGradient, lane B
// 2026-11-13), against the central difference of the CACHED kernels it differentiates (GaussianRF::Overlap2C / Grad2 / Repulsion2C with the ket radial
// placed at B ± h), per shell pair and weight set, on every Cartesian shell pair through d and on contracted radials.  This is the gate
// at the kernel itself: the lattice pullback (M_LatticePullback) sums it over images and atoms, where a wrong per-image derivative can
// cancel in a gradient and show only in the virial.
#include "gtest/gtest.h"
#include <cmath>
#include <cstdio>
#include <vector>

import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;     // GaussianRF (a test may cheat: the kernel itself)
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Structure;                                               // Molecule, Site (the nuclei)
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;

namespace {

std::vector<Polarization> ShellComponents(int L)
{
    std::vector<Polarization> p;
    for (int n=L;n>=0;n--) for (int l=L-n;l>=0;l--) p.push_back(Polarization(n,l,L-n-l));
    return p;
}

//! f(B) = Σ_pq W_pq K(ra, pa[p]; rb@B, pb[q]) by the CACHED named kernels.
double Value(GaussianRF::TwoCentre K, const GaussianRF& ra, const GaussianRF& rb, const rvec3_t& B,
             const std::vector<Polarization>& pa, const std::vector<Polarization>& pb, const std::vector<double>& W)
{
    const GaussianRF rbB=rb.AtCenter(B);
    double s=0.0;
    for (size_t p=0;p<pa.size();p++)
        for (size_t q=0;q<pb.size();q++)
            s+=W[p*pb.size()+q]*(K==GaussianRF::TwoCentre::Overlap ? ra.Overlap2C(rbB,pa[p],pb[q])
                                : K==GaussianRF::TwoCentre::Grad2   ? ra.Grad2(rbB,pa[p],pb[q]) : ra.Repulsion2C(rbB,pa[p],pb[q]));
    return s;
}

double Comp(const rvec3_t& v, int c) {return c==0 ? v.x : c==1 ? v.y : v.z;}
const char* KernelName(GaussianRF::TwoCentre K) {return K==GaussianRF::TwoCentre::Overlap ? "S" : K==GaussianRF::TwoCentre::Grad2 ? "T" : "(a|b)";}

} // namespace

TEST(KetGradient, IsTheFiniteDifferenceOfTheCachedKernelForEveryShellPairThroughD)
{
    const rvec3_t A(0.13,-0.21,0.34);
    double worst=0.0;
    // A generic ket centre, and two where values vanish BY SYMMETRY (on an axis through the bra, and on a coordinate plane through it):
    // a kernel that branches on a zero value has a zero derivative there that is not zero (the 2026-11-13 lattice-virial failure).  And
    // the ket ON the bra (B = A): the Coulomb metric's same-atom pairs, where T = 0 in the Boys function.
    for (const rvec3_t& B0 : {rvec3_t(1.07,0.52,-0.61), rvec3_t(A.x+1.3,A.y,A.z), rvec3_t(A.x+0.9,A.y-0.8,A.z), A})
    for (int contracted=0;contracted<2;contracted++)
        for (int La=0;La<=2;La++)
            for (int Lb=0;Lb<=2;Lb++)
            {
                const GaussianRF ra = contracted ? GaussianRF(rvec_t{0.6,0.4}, rvec_t{1.9,0.45}, A, La) : GaussianRF(0.83, A, La);
                const GaussianRF rb = contracted ? GaussianRF(rvec_t{0.3,0.7}, rvec_t{2.6,0.38}, B0, Lb) : GaussianRF(0.57, B0, Lb);
                const std::vector<Polarization> pa=ShellComponents(La), pb=ShellComponents(Lb);
                std::vector<double> W(pa.size()*pb.size());
                for (size_t k=0;k<W.size();k++) W[k]=std::sin(1.3*double(k)+0.7*La+0.3*Lb+contracted);   // dense, signed, no zeros
                for (auto K : {GaussianRF::TwoCentre::Overlap, GaussianRF::TwoCentre::Grad2, GaussianRF::TwoCentre::Repulsion})
                {
                    const rvec3_t g=ra.KetGradient(K, rb, B0, pa.data(), pa.size(), pb.data(), pb.size(), W.data());
                    for (int c=0;c<3;c++)
                    {
                        const double h=1e-4;
                        auto at=[&](double s){ rvec3_t B=B0; (c==0 ? B.x : c==1 ? B.y : B.z)+=s; return Value(K,ra,rb,B,pa,pb,W); };
                        const double fd=(4.0*(at(0.5*h)-at(-0.5*h))/h - (at(h)-at(-h))/(2*h))/3.0;   // Richardson on the central difference
                        const double an=Comp(g,c), err=std::fabs(an-fd), tol=1e-8*std::max(1.0,std::fabs(fd));
                        worst=std::max(worst,err);
                        EXPECT_LT(err, tol) << KernelName(K) << " La=" << La << " Lb=" << Lb
                                            << " contracted=" << contracted << " B-A=(" << B0.x-A.x << "," << B0.y-A.y << "," << B0.z-A.z << ") axis " << c << ": AD " << an << " FD " << fd;
                    }
                }
            }
    std::printf("  [ket gradient] worst |AD - FD| over S, <p^2> and (a|b), La,Lb <= 2, contracted and not, generic, on-axis, on-plane and coincident ket: %.2e\n", worst);
}

// THE NUCLEAR TWIN (GaussianRF::NuclearGradient), every centre moving, against the central difference of the CACHED Nuclear: the bra, the
// ket and each nucleus displaced in turn.  Two nuclei, one on the bra-ket axis (values vanishing by symmetry: the trap a zero-value skip
// falls into) and one off it; every Cartesian shell pair through d, contracted and not.
TEST(KetGradient, TheNuclearTwinIsTheFiniteDifferenceOfTheCachedKernelWithEveryCentreMoving)
{
    const rvec3_t A0(0.13,-0.21,0.34), B0(A0.x+1.3,A0.y,A0.z);
    const std::vector<rvec3_t> C0={rvec3_t(A0.x+0.6,A0.y,A0.z), rvec3_t(0.9,0.8,-0.5)};
    const std::vector<double>  Z ={3.0, 1.0};
    auto mol=[&](const std::vector<rvec3_t>& C){ Molecule m; for (size_t k=0;k<C.size();k++) m.Insert(Site{int(Z[k]),C[k]}); return m; };
    double worst=0.0;
    for (int contracted=0;contracted<2;contracted++)
        for (int La=0;La<=2;La++)
            for (int Lb=0;Lb<=2;Lb++)
            {
                auto radA=[&](const rvec3_t& A){ return contracted ? GaussianRF(rvec_t{0.6,0.4}, rvec_t{1.9,0.45}, A, La) : GaussianRF(0.83, A, La); };
                auto radB=[&](const rvec3_t& B){ return contracted ? GaussianRF(rvec_t{0.3,0.7}, rvec_t{2.6,0.38}, B, Lb) : GaussianRF(0.57, B, Lb); };
                const std::vector<Polarization> pa=ShellComponents(La), pb=ShellComponents(Lb);
                std::vector<double> W(pa.size()*pb.size());
                for (size_t k=0;k<W.size();k++) W[k]=std::cos(0.9*double(k)+0.4*La+0.2*Lb+contracted);
                auto value=[&](const rvec3_t& A, const rvec3_t& B, const std::vector<rvec3_t>& C)
                {
                    const GaussianRF ra=radA(A), rb=radB(B); const Molecule m=mol(C);
                    double v=0.0;
                    for (size_t p=0;p<pa.size();p++) for (size_t q=0;q<pb.size();q++) v+=W[p*pb.size()+q]*ra.Nuclear(rb,pa[p],pb[q],&m);
                    return v;
                };
                rvec3_t gA, gB; std::vector<rvec3_t> gC;
                radA(A0).NuclearGradient(radB(B0), pa.data(), pa.size(), pb.data(), pb.size(), W.data(), C0, Z, gA, gB, gC);
                ASSERT_EQ(gC.size(), C0.size());
                // centre 0 = A, 1 = B, 2.. = the nuclei
                for (size_t centre=0;centre<2+C0.size();centre++)
                    for (int c=0;c<3;c++)
                    {
                        const double h=1e-4;
                        auto at=[&](double s)
                        {
                            rvec3_t A=A0, B=B0; std::vector<rvec3_t> C=C0;
                            rvec3_t& X = centre==0 ? A : centre==1 ? B : C[centre-2];
                            (c==0 ? X.x : c==1 ? X.y : X.z)+=s;
                            return value(A,B,C);
                        };
                        const double fd=(4.0*(at(0.5*h)-at(-0.5*h))/h - (at(h)-at(-h))/(2*h))/3.0;
                        const rvec3_t& g = centre==0 ? gA : centre==1 ? gB : gC[centre-2];
                        const double an=Comp(g,c), err=std::fabs(an-fd), tol=1e-8*std::max(1.0,std::fabs(fd));
                        worst=std::max(worst,err);
                        EXPECT_LT(err, tol) << "La=" << La << " Lb=" << Lb << " contracted=" << contracted << " centre " << centre
                                            << " axis " << c << ": AD " << an << " FD " << fd;
                    }
            }
    std::printf("  [nuclear twin] worst |AD - FD| over bra, ket and two nuclei (one on the axis), La,Lb <= 2: %.2e\n", worst);
}


// THE 3-CENTRE TWINS (GaussianRF::ThreeCentreGradient), every centre moving, against the central difference of the CACHED Overlap3C and
// Repulsion3C (this radial = centre C, as they have it).  C on the A-B axis (Hermite coefficients and RNLM entries vanishing by symmetry:
// the trap the cached kernels' zero-value skips fall into) and off it; every Cartesian pair through d; the field radial of L 0 and 2 with
// two terms (a polynomial field, as a GTH local part's r^{2n} terms make).
TEST(KetGradient, TheThreeCentreTwinsAreTheFiniteDifferenceOfTheCachedKernelsWithEveryCentreMoving)
{
    const rvec3_t A0(0.13,-0.21,0.34), B0(A0.x+1.3,A0.y,A0.z);
    double worst[2]={0.0,0.0};
    for (const rvec3_t& C0 : {rvec3_t(A0.x+0.5,A0.y,A0.z), rvec3_t(0.9,0.8,-0.5)})
    for (int Lc : {0,2})
    for (int La=0;La<=2;La++)
        for (int Lb=0;Lb<=2;Lb++)
        {
            auto radA=[&](const rvec3_t& A){ return GaussianRF(rvec_t{0.6,0.4}, rvec_t{1.9,0.45}, A, La); };
            auto radB=[&](const rvec3_t& B){ return GaussianRF(0.57, B, Lb); };
            auto radC=[&](const rvec3_t& C){ return GaussianRF(0.71, C, Lc); };
            const std::vector<Polarization> pa=ShellComponents(La), pb=ShellComponents(Lb);
            const std::vector<Polarization> pt = Lc==0 ? std::vector<Polarization>{Polarization(0,0,0)} : std::vector<Polarization>{Polarization(2,0,0), Polarization(0,1,1)};
            const std::vector<double> ct = Lc==0 ? std::vector<double>{1.0} : std::vector<double>{0.8,-0.3};
            std::vector<double> W(pa.size()*pb.size());
            for (size_t k=0;k<W.size();k++) W[k]=std::sin(0.7*double(k)+0.5*La+0.25*Lb+Lc);
            for (auto K : {GaussianRF::ThreeCentre::Overlap, GaussianRF::ThreeCentre::Repulsion})
            {
                auto value=[&](const rvec3_t& A, const rvec3_t& B, const rvec3_t& C)
                {
                    const GaussianRF ra=radA(A), rb=radB(B), rc=radC(C);
                    double v=0.0;
                    for (size_t p=0;p<pa.size();p++) for (size_t q=0;q<pb.size();q++) for (size_t u=0;u<pt.size();u++)
                        v+=W[p*pb.size()+q]*ct[u]*(K==GaussianRF::ThreeCentre::Overlap ? rc.Overlap3C(ra,rb,pa[p],pb[q],pt[u])
                                                                                        : rc.Repulsion3C(ra,rb,pa[p],pb[q],pt[u]));
                    return v;
                };
                rvec3_t gA, gB, gC;
                radC(C0).ThreeCentreGradient(K, radA(A0), radB(B0), pa.data(), pa.size(), pb.data(), pb.size(), pt.data(), ct.data(), pt.size(), W.data(), gA, gB, gC);
                for (int centre=0;centre<3;centre++)
                    for (int c=0;c<3;c++)
                    {
                        const double h=1e-4;
                        auto at=[&](double s)
                        {
                            rvec3_t A=A0, B=B0, C=C0;
                            rvec3_t& X = centre==0 ? A : centre==1 ? B : C;
                            (c==0 ? X.x : c==1 ? X.y : X.z)+=s;
                            return value(A,B,C);
                        };
                        const double fd=(4.0*(at(0.5*h)-at(-0.5*h))/h - (at(h)-at(-h))/(2*h))/3.0;
                        const rvec3_t& g = centre==0 ? gA : centre==1 ? gB : gC;
                        const double an=Comp(g,c), err=std::fabs(an-fd), tol=1e-8*std::max(1.0,std::fabs(fd));
                        double& wst=worst[K==GaussianRF::ThreeCentre::Overlap ? 0 : 1];
                        wst=std::max(wst,err);
                        EXPECT_LT(err, tol) << (K==GaussianRF::ThreeCentre::Overlap ? "<ab|c>" : "(ab|c)") << " La=" << La << " Lb=" << Lb << " Lc=" << Lc
                                            << " centre " << centre << " axis " << c << ": AD " << an << " FD " << fd;
                    }
            }
        }
    std::printf("  [3-centre twins] worst |AD - FD|: <ab|c> %.2e, (ab|c) %.2e\n", worst[0], worst[1]);
}

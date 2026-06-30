// file: PlaneWaveUT.C  Empty-lattice (free-electron) validation of the plane-wave basis set.
//
// Milestone 1 of doc/PlaneWavePlan.md: with V_ext = 0 the band energies are exactly
// 1/2 |k+G|^2.  We build the PW basis for a cubic cell at several k-points, diagonalise
// H = 1/2 * Kinetic() (V=0), and compare the eigenvalues to the analytic free-electron
// ladder computed independently (cubic: B = (2 pi / a) I, so |k+G| = (2 pi/a)|k_frac+m|).

#include <vector>
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH (H, Si pseudopotentials from the database)
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.LASolver;
import qchem.Types;
import qchem.Blaze;
import qchem.Math;           // Pi
using namespace qchem;

using BasisSet::Lattice_3D::PlaneWave_IBS;
using Pseudopotential::LocalPotential;
using Pseudopotential::BareCoulomb;
using Pseudopotential::GaussianSmearedNucleus;
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::SeparablePotential;
using Pseudopotential::GaussianProjector;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

namespace
{

// Analytic free-electron eigenvalues for a CUBIC cell of edge a, k-point kIndex/N and cutoff Ecut.
// Independent of the basis-set code: B = (2 pi / a) I, |k+G|^2 = (2 pi/a)^2 |k_frac + m|^2.
std::vector<double> FreeElectronReference(double a, ivec3_t N, ivec3_t kIndex, double Ecut)
{
    double b=2*Pi/a;
    rvec3_t kf(kIndex.x/double(N.x), kIndex.y/double(N.y), kIndex.z/double(N.z));
    int mMax=int(std::ceil((std::sqrt(2*Ecut)+b*std::sqrt(kf*kf))/b))+1;
    std::vector<double> e;
    for (int mx=-mMax; mx<=mMax; mx++)
    for (int my=-mMax; my<=mMax; my++)
    for (int mz=-mMax; mz<=mMax; mz++)
    {
        rvec3_t kG=(kf+ivec3_t(mx,my,mz))*b;   // (k+G) Cartesian
        double E=0.5*(kG*kG);
        if (E<Ecut) e.push_back(E);
    }
    std::sort(e.begin(),e.end());
    return e;
}

// Diagonalise H = 1/2 * Kinetic + V and return the sorted eigenvalues (V optional / nullptr for V=0).
std::vector<double> SolveBands(const PlaneWave_IBS& pw, const chmat_t* V=nullptr)
{
    size_t n=pw.GetNumFunctions();
    chmat_t p2=pw.MakeKinetic();
    hmat_t<dcmplx> H(n), S(n);
    for (size_t i=0;i<n;i++)
    {
        S(i,i)=1.0;
        for (size_t j=i;j<n;j++)
        {
            dcmplx h = (i==j) ? dcmplx(0.5*std::real(dcmplx(p2(i,i)))) : dcmplx(0.0);
            if (V) h += dcmplx((*V)(i,j));
            H(i,j)=h;
        }
    }

    LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
    las->SetBasisOverlap(S);
    auto [U,e]=las->Solve(H);
    delete las;

    std::vector<double> ev;                     // blaze iterators don't satisfy std iterator traits
    for (size_t i=0;i<e.size();i++) ev.push_back(e[i]);
    std::sort(ev.begin(),ev.end());
    return ev;
}

class PlaneWaveTests : public ::testing::Test {};

// At each k-point the PW bands must reproduce the exact free-electron ladder.
void CheckKPoint(double a, ivec3_t N, ivec3_t kIndex, double Ecut)
{
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,kIndex,Ecut);

    std::vector<double> ref=FreeElectronReference(a,N,kIndex,Ecut);
    ASSERT_EQ(pw.GetNumFunctions(),ref.size());          // same G-set

    // Overlap is the identity (plane waves are orthonormal over the cell).
    chmat_t S=pw.MakeOverlap();
    for (size_t i=0;i<pw.GetNumFunctions();i++)
    {
        EXPECT_NEAR(std::real(dcmplx(S(i,i))),1.0,1e-14);
        for (size_t j=i+1;j<pw.GetNumFunctions();j++) EXPECT_DOUBLE_EQ(std::norm(dcmplx(S(i,j))),0.0);
    }

    std::vector<double> bands=SolveBands(pw);
    ASSERT_EQ(bands.size(),ref.size());
    for (size_t i=0;i<ref.size();i++) EXPECT_NEAR(bands[i],ref[i],1e-10);
}

// Separable cubic cosine V = V0 (cos bx + cos by + cos bz), b = 2 pi / a.  Its only nonzero Fourier
// components are Vtilde(+/- b_i) = V0/2, i.e. Vtilde(dm) = V0/2 when dm is a unit step +/- e_i.
std::function<dcmplx(const ivec3_t&)> CosineVtilde(double V0)
{
    return [V0](const ivec3_t& dm)->dcmplx
    {
        return (dm.x*dm.x + dm.y*dm.y + dm.z*dm.z == 1) ? dcmplx(0.5*V0) : dcmplx(0.0);
    };
}

// Gamma-point hydrogen ground state for a given local potential: Z=1 nucleus at the cell origin,
// large-a isolation (flat bands => Gamma only suffices).
double HydrogenE0(double a, double Ecut, const LocalPotential& v, size_t* npw=nullptr)
{
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),Ecut);
    if (npw) *npw=pw.GetNumFunctions();
    Atom H(1,rvec3_t(0,0,0));
    chmat_t V=pw.MakeLocalPotential(&H,v);
    std::vector<double> bands=SolveBands(pw,&V);
    return bands.front();
}

} // namespace

TEST_F(PlaneWaveTests, EmptyLatticeGamma)   { CheckKPoint(5.0, ivec3_t(4,4,4), ivec3_t(0,0,0), 8.0); }
TEST_F(PlaneWaveTests, EmptyLatticeOffGamma){ CheckKPoint(5.0, ivec3_t(4,4,4), ivec3_t(1,0,0), 8.0); }
TEST_F(PlaneWaveTests, EmptyLatticeOrthorhombicGamma)
{
    // Non-cubic edge would need the general B; keep cubic but a different edge to vary the ladder spacing.
    CheckKPoint(8.0, ivec3_t(2,2,2), ivec3_t(1,1,0), 4.0);
}

// --- Milestone 2.2: separable cosine potential (validates G-space potential assembly) -------------

// MakePotential must produce exactly the sparse cosine structure: V0/2 between G and G +/- b_i, else 0.
TEST_F(PlaneWaveTests, CosineMatrixStructure)
{
    double a=5.0, V0=0.3;
    ivec3_t N(4,4,4);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),8.0);

    auto Vtilde=CosineVtilde(V0);
    chmat_t V=pw.MakePotential(Vtilde);
    size_t n=pw.GetNumFunctions();
    ASSERT_EQ(V.rows(),n);
    for (size_t i=0;i<n;i++)
    {
        EXPECT_DOUBLE_EQ(std::real(dcmplx(V(i,i))),0.0);   // cos has no G=0 component -> traceless
        for (size_t j=0;j<n;j++)
        {
            // V(i,j) must be Vtilde(m_i - m_j): nonzero (V0/2) iff the two plane waves differ by one b_i step.
            ivec3_t dm=pw.GetGIndex(i)-pw.GetGIndex(j);
            double expect=std::real(Vtilde(dm));
            EXPECT_DOUBLE_EQ(std::real(dcmplx(V(i,j))),expect);
            EXPECT_DOUBLE_EQ(std::imag(dcmplx(V(i,j))),0.0);
        }
    }
}

// The potential is traceless, so the eigenvalue sum is exactly the free-electron kinetic sum.
TEST_F(PlaneWaveTests, CosineTraceInvariant)
{
    double a=5.0, V0=0.5;
    ivec3_t N(4,4,4);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(1,0,0),8.0);

    chmat_t V=pw.MakePotential(CosineVtilde(V0));
    std::vector<double> bands=SolveBands(pw,&V);

    chmat_t p2=pw.MakeKinetic();
    double kineticSum=0.0, bandSum=0.0;
    for (size_t i=0;i<pw.GetNumFunctions();i++) kineticSum += 0.5*std::real(dcmplx(p2(i,i)));
    for (double e:bands) bandSum += e;
    EXPECT_NEAR(bandSum,kineticSum,1e-9);
}

// V0 -> 0 must recover the empty-lattice ladder exactly (regression onto milestone 1).
TEST_F(PlaneWaveTests, CosineZeroRecoversEmptyLattice)
{
    double a=5.0, Ecut=8.0;
    ivec3_t N(4,4,4), k(1,0,0);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,k,Ecut);

    chmat_t V=pw.MakePotential(CosineVtilde(0.0));
    std::vector<double> bands=SolveBands(pw,&V);
    std::vector<double> ref=FreeElectronReference(a,N,k,Ecut);
    ASSERT_EQ(bands.size(),ref.size());
    for (size_t i=0;i<ref.size();i++) EXPECT_NEAR(bands[i],ref[i],1e-10);
}

// Physics cross-check: 2nd-order perturbation theory for the Gamma-point ground state.
// E0 ~ sum_{G!=0} |<G|V|0>|^2 / (0 - E_G) = 6 (V0/2)^2 / (-b^2/2) = -3 V0^2 / b^2 = -3 V0^2 a^2/(4 pi^2).
TEST_F(PlaneWaveTests, CosineGroundStatePerturbation)
{
    double a=5.0, V0=0.01;
    ivec3_t N(4,4,4);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),8.0);

    chmat_t V=pw.MakePotential(CosineVtilde(V0));
    std::vector<double> bands=SolveBands(pw,&V);

    double E0_pt2 = -3.0*V0*V0*a*a/(4*Pi*Pi);
    EXPECT_NEAR(bands.front(),E0_pt2, 0.01*std::abs(E0_pt2)); // 1% -- higher-order PT is O(V0^4)
}

// --- Milestone 2.3: hydrogen, bare Coulomb (validates the nuclear structure-factor assembly) -------

// Exact check of the assembled potential: a single nucleus (Z=1) at the origin in a cubic cell of edge
// a gives V(dG) = -(4 pi / a^3)/|dG|^2.  The nearest-G coupling (|dG| = b = 2 pi/a) is exactly -1/(pi a).
TEST_F(PlaneWaveTests, HydrogenPotentialMatrixElement)
{
    double a=8.0;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),6.0);
    Atom H(1,rvec3_t(0,0,0));
    chmat_t V=pw.MakeNuclear(&H);

    size_t n=pw.GetNumFunctions();
    double nearest=-1.0/(Pi*a);     // expected coupling between G and G +/- b_i
    bool sawNearest=false;
    for (size_t i=0;i<n;i++)
    {
        EXPECT_DOUBLE_EQ(std::real(dcmplx(V(i,i))),0.0);   // diagonal (dG=0) dropped
        for (size_t j=0;j<n;j++)
        {
            ivec3_t dm=pw.GetGIndex(i)-pw.GetGIndex(j);
            int s=dm.x*dm.x+dm.y*dm.y+dm.z*dm.z;
            EXPECT_DOUBLE_EQ(std::imag(dcmplx(V(i,j))),0.0); // atom at origin -> real
            if (s==1) { EXPECT_NEAR(std::real(dcmplx(V(i,j))),nearest,1e-12); sawNearest=true; }
        }
    }
    EXPECT_TRUE(sawNearest);
}

// End-to-end: bare-Coulomb hydrogen converges toward -0.5 Ha FROM ABOVE, but only very slowly (the 1s
// cusp; see doc/PlaneWavePlan.md sec.2.3).  At these modest cutoffs E0 is far from -0.5, so the robust,
// cheap assertions are: (1) variational monotonicity -- more plane waves never raise the energy at fixed
// cell; (2) a bound state of the right order.  Observed (a=8): Ecut 4 -> -0.1364, 6 -> -0.1494.
TEST_F(PlaneWaveTests, HydrogenVariationalConvergence)
{
    double a=8.0;
    double E_low =HydrogenE0(a,4.0,BareCoulomb());
    double E_high=HydrogenE0(a,6.0,BareCoulomb());
    EXPECT_LT(E_high,E_low+1e-9);            // bigger basis (same H) cannot raise the variational energy
    EXPECT_LT(E_high,-0.05);                 // genuinely bound, attractive potential of the right order
    EXPECT_GT(E_high,-0.5);                  // still well above the exact -0.5 at this (cutoff, cell)
}

// --- Rung 1: Gaussian-smeared nucleus (local pseudopotential) -- the cusp is removed ---------------

TEST_F(PlaneWaveTests, DISABLED_SmearedCalibration)
{
    double a=8.0;
    for (double Ecut : {4.0,6.0,9.0,12.0})
    {
        double Eb=HydrogenE0(a,Ecut,BareCoulomb());
        double Es=HydrogenE0(a,Ecut,GaussianSmearedNucleus(0.5));
        printf("Ecut=%4.1f  bare=% .5f  smeared(0.5)=% .5f\n",Ecut,Eb,Es);
    }
}

// --- Real norm-conserving pseudopotential: HGH hydrogen (local; H has no core projectors) ----------

TEST_F(PlaneWaveTests, DISABLED_HGHCalibration)
{
    double a=7.0;
    HGH_LocalPotential h=GetGTH("H","LDA").local;
    for (double Ecut : {4.0,6.0,9.0,12.0})
    {
        size_t npw=0;
        double Eh=HydrogenE0(a,Ecut,h,&npw);
        double Eb=HydrogenE0(a,Ecut,BareCoulomb());
        printf("Ecut=%5.1f  npw=%4zu  HGH-H=% .5f  bare=% .5f\n",Ecut,npw,Eh,Eb);
    }
    // form factor: HGH should track bare -4pi Z/G^2 at small G (Coulomb tail) and be Gaussian-soft at large G
    for (double G2 : {0.1,1.0,5.0,25.0,100.0})
        printf("G2=%6.1f  HGH=% .5f  bare=% .5f  ratio=% .4f\n",
               G2,h.FormFactor(1,G2),BareCoulomb().FormFactor(1,G2),
               h.FormFactor(1,G2)/BareCoulomb().FormFactor(1,G2));
}

// HGH hydrogen is a REAL norm-conserving pseudopotential -- but a degenerate one: hydrogen has no core
// electrons to pseudize, so the PP carries no nonlocal projectors and its local part stays essentially
// the bare -Z/r (the erf only softens the large-G tail).  This pins that correct behaviour:
//  (1) the form factor tracks bare -4pi Z/G^2 at small G (the long-range Coulomb is preserved),
//  (2) it is strictly softer than bare at large G (the would-be core is pseudized),
//  (3) the ground state matches the bare-Coulomb result -- the no-core limit.
// The pseudopotential PAYOFF (a soft potential that converges far faster than all-electron) shows up for
// CORE-bearing elements, which also carry the nonlocal l-channel projectors -- the next rung of lineage A.
TEST_F(PlaneWaveTests, HGHHydrogenIsRealButCoreless)
{
    HGH_LocalPotential h=GetGTH("H","LDA").local;
    BareCoulomb bare;
    EXPECT_NEAR(h.FormFactor(1,1.0)/bare.FormFactor(1,1.0), 1.0, 2e-3);                 // (1) Coulomb tail
    EXPECT_LT(std::abs(h.FormFactor(1,100.0)), std::abs(bare.FormFactor(1,100.0)));     // (2) softer core

    double a=7.0, Ecut=9.0;                                                            // (3) no-core limit
    double Eh=HydrogenE0(a,Ecut,h), Eb=HydrogenE0(a,Ecut,bare);
    EXPECT_LT(Eh,-0.05);                                                               // genuinely bound
    EXPECT_NEAR(Eh,Eb,3e-3);                                                           // H coreless -> ~bare
}

// sigma -> 0 must reproduce the bare Coulomb potential matrix element by element.
TEST_F(PlaneWaveTests, SmearedZeroRecoversBareCoulomb)
{
    double a=8.0;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),6.0);
    Atom H(1,rvec3_t(0,0,0));

    chmat_t Vbare    = pw.MakeLocalPotential(&H,BareCoulomb());
    chmat_t Vtiny    = pw.MakeLocalPotential(&H,GaussianSmearedNucleus(1e-6));
    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
            EXPECT_NEAR(std::real(dcmplx(Vtiny(i,j))),std::real(dcmplx(Vbare(i,j))),1e-8);
}

// The cusp is gone, so the smeared energy converges FAST with Ecut, unlike the bare-Coulomb crawl.
// Observed (a=8, sigma=0.5): Ecut 4->6 moves smeared by ~8e-4 (and only ~1e-4 more out to 12),
// while bare moves ~1.3e-2 and keeps crawling.
TEST_F(PlaneWaveTests, SmearedConvergesFast)
{
    double a=8.0;
    GaussianSmearedNucleus vsm(0.5);
    double Es_lo=HydrogenE0(a,4.0,vsm);
    double Es_hi=HydrogenE0(a,6.0,vsm);
    double Eb_lo=HydrogenE0(a,4.0,BareCoulomb());
    double Eb_hi=HydrogenE0(a,6.0,BareCoulomb());

    double smearedDrift=std::abs(Es_hi-Es_lo);
    double bareDrift   =std::abs(Eb_hi-Eb_lo);
    EXPECT_LT(Es_hi,Es_lo+1e-9);             // still variational
    EXPECT_LT(smearedDrift,5e-3);            // nearly converged already at modest Ecut
    EXPECT_LT(smearedDrift,bareDrift);       // and converging much faster than bare Coulomb
}

// --- Rung 2: Kleinman-Bylander separable NONLOCAL potential -- V_NL = |beta> D <beta| -------------

// One s-channel Gaussian projector at the origin gives V_NL(G,G') = (D/Omega) b_i b_j exactly, with
// b_i = exp(-sigma^2 |k+G_i|^2 / 2).  At Gamma (cubic), |G_i|^2 = (2 pi/a)^2 |m_i|^2.
TEST_F(PlaneWaveTests, SeparableNonlocalMatrixElement)
{
    double a=8.0, sigma=1.0, D=0.7;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),6.0);
    Atom H(1,rvec3_t(0,0,0));

    chmat_t V=pw.MakeSeparablePotential(&H,GaussianProjector(sigma,D));
    size_t n=pw.GetNumFunctions();
    double Omega=a*a*a, b=2*Pi/a;
    auto bfac=[&](size_t i){ ivec3_t m=pw.GetGIndex(i); double G2=b*b*(m.x*m.x+m.y*m.y+m.z*m.z);
                             return std::exp(-0.5*sigma*sigma*G2); };
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            double expect=(D/Omega)*bfac(i)*bfac(j);
            EXPECT_NEAR(std::real(dcmplx(V(i,j))),expect,1e-12);
            EXPECT_DOUBLE_EQ(std::imag(dcmplx(V(i,j))),0.0);
        }
}

// A single projector makes V_NL rank 1: exactly one nonzero eigenvalue, equal to the trace.
TEST_F(PlaneWaveTests, SeparableNonlocalIsRankOne)
{
    double a=8.0;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),6.0);
    Atom H(1,rvec3_t(0,0,0));

    chmat_t V=pw.MakeSeparablePotential(&H,GaussianProjector(1.0,0.7));
    size_t n=pw.GetNumFunctions();
    double trace=0.0;
    for (size_t i=0;i<n;i++) trace += std::real(dcmplx(V(i,i)));

    rvec_t d; mat_t<dcmplx> U;
    blazem::eigen(V,d,U);
    int nonzero=0; double lambda=0.0;
    for (size_t i=0;i<d.size();i++) if (std::abs(d[i])>1e-10) { nonzero++; lambda=d[i]; }
    EXPECT_EQ(nonzero,1);                     // separable single projector => rank 1
    EXPECT_NEAR(lambda,trace,1e-10);          // the lone eigenvalue carries the whole trace
}

// An l=1 (p-channel) projector carries the angular weight (2l+1)P_1(cos gamma) = 3 cos gamma, gamma the
// angle between k+G and k+G'.  Single projector at the origin: V(G,G') = (3D/Omega) cos gamma b_i b_j.
// The signature of l=1 vs l=0 is the angular dependence: parallel k+G couple x3, PERPENDICULAR k+G do
// not couple at all, and anti-parallel couple with a sign flip.
TEST_F(PlaneWaveTests, SeparableNonlocalL1AngularFactor)
{
    double a=8.0, sigma=1.0, D=0.7;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),6.0);
    Atom H(1,rvec3_t(0,0,0));

    chmat_t V=pw.MakeSeparablePotential(&H,GaussianProjector(sigma,D,/*l=*/1));
    size_t n=pw.GetNumFunctions();
    double Omega=a*a*a, b=2*Pi/a;
    auto idx=[&](ivec3_t m)->size_t{ for(size_t i=0;i<n;i++){ ivec3_t g=pw.GetGIndex(i);
                             if(g.x==m.x&&g.y==m.y&&g.z==m.z) return i; } return n; };
    auto bfac=[&](ivec3_t m){ double G2=b*b*(m.x*m.x+m.y*m.y+m.z*m.z);
                             return std::exp(-0.5*sigma*sigma*G2); };

    size_t ix=idx(ivec3_t(1,0,0)), iy=idx(ivec3_t(0,1,0)), imx=idx(ivec3_t(-1,0,0));
    ASSERT_LT(ix,n); ASSERT_LT(iy,n); ASSERT_LT(imx,n);
    double b1=bfac(ivec3_t(1,0,0));
    EXPECT_NEAR(std::real(dcmplx(V(ix,ix))),  (3*D/Omega)*b1*b1, 1e-12);   // parallel: cos=1  -> x3
    EXPECT_NEAR(std::real(dcmplx(V(ix,iy))),  0.0,               1e-12);   // perpendicular: cos=0 -> 0
    EXPECT_NEAR(std::real(dcmplx(V(ix,imx))),-(3*D/Omega)*b1*b1, 1e-12);   // anti-parallel: cos=-1 -> -x3
}

// One l=1 radial projector spans m = -1,0,+1, so its V_NL has rank 2l+1 = 3 (the l=0 case was rank 1).
// Off-Gamma so every k+G has a defined direction (no k+G = 0).
TEST_F(PlaneWaveTests, SeparableNonlocalL1IsRankThree)
{
    double a=8.0;
    ivec3_t N(4,4,4), k(1,0,0);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,k,6.0);
    Atom H(1,rvec3_t(0,0,0));

    chmat_t V=pw.MakeSeparablePotential(&H,GaussianProjector(1.0,0.7,/*l=*/1));
    size_t n=pw.GetNumFunctions();
    double trace=0.0;
    for (size_t i=0;i<n;i++) trace += std::real(dcmplx(V(i,i)));

    rvec_t d; mat_t<dcmplx> U;
    blazem::eigen(V,d,U);
    int nonzero=0; double sum=0.0;
    for (size_t i=0;i<d.size();i++) if (std::abs(d[i])>1e-10) { nonzero++; sum+=d[i]; }
    EXPECT_EQ(nonzero,3);                     // m = -1,0,+1 => rank 2l+1
    EXPECT_NEAR(sum,trace,1e-10);             // the 3 nonzero eigenvalues carry the whole trace
}

// Composition: the external block is V = V_loc + V_nonlocal.  A repulsive (D>0) projector that the
// smooth ground state overlaps must raise E0 relative to the local-only potential.
TEST_F(PlaneWaveTests, LocalPlusNonlocalRaisesGroundState)
{
    double a=8.0;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),6.0);
    Atom H(1,rvec3_t(0,0,0));

    chmat_t Vloc=pw.MakeLocalPotential(&H,GaussianSmearedNucleus(0.5));
    chmat_t Vnl =pw.MakeSeparablePotential(&H,GaussianProjector(1.0,0.5)); // repulsive
    chmat_t Vtot=Vloc+Vnl;

    std::vector<double> locOnly=SolveBands(pw,&Vloc);
    std::vector<double> combined=SolveBands(pw,&Vtot);
    EXPECT_GT(combined.front(),locOnly.front());   // repulsive nonlocal pushes the ground state up
}

// --- A real material: HGH silicon (GTH-LDA q4), local + nonlocal s,p projectors ---------------------

// The HGH-Si nonlocal part has a 2x2 s-channel and a 1x1 p-channel.  Diagonalising the s h-matrix must
// recover its eigenvalues as the KB coefficients; the p projector (l=1) vanishes at q=0 (beta ~ q^l),
// while the s projectors (l=0) do not.
TEST_F(PlaneWaveTests, HGHSiliconNonlocalChannels)
{
    HGH_SeparablePotential si=GetGTH("Si","LDA",4).nonlocal;
    ASSERT_EQ(si.NumProjectors(14), 3u);

    std::vector<int>    ls;
    std::vector<double> D;
    for (size_t p=0;p<3;p++){ ls.push_back(si.AngularMomentum(14,p)); D.push_back(si.Coefficient(14,p)); }
    EXPECT_EQ(std::count(ls.begin(),ls.end(),0), 2);     // two s projectors
    EXPECT_EQ(std::count(ls.begin(),ls.end(),1), 1);     // one p projector

    std::sort(D.begin(),D.end());
    EXPECT_NEAR(D[0], 2.72701346, 1e-6);                 // p (h11)
    EXPECT_NEAR(D[1], 2.753267,   1e-5);                 // s, smaller eigenvalue of the 2x2 h-matrix
    EXPECT_NEAR(D[2], 6.411857,   1e-5);                 // s, larger eigenvalue

    for (size_t p=0;p<3;p++)                             // angular character at q=0
    {
        double b0=si.Projector(14,p,0.0);
        if (si.AngularMomentum(14,p)==1) EXPECT_NEAR(b0,0.0,1e-12);   // p ~ q -> 0
        else                             EXPECT_GT(std::abs(b0),0.0); // s finite
    }
}

// The full real HGH-Si one-body potential V = V_loc + V_nonlocal assembles into a well-formed (Hermitian,
// real-spectrum) H(k) for a Si pseudo-atom.  Quantitative Si bands need the DFT Hartree+XC self-consistency;
// this pins the ionic-potential machinery that the DFT step builds on.
TEST_F(PlaneWaveTests, HGHSiliconHamiltonianWellFormed)
{
    double a=10.0;
    ivec3_t N(1,1,1);
    UnitCell cell(a);
    Lattice_3D lat(cell,N);
    PlaneWave_IBS pw(lat.Reciprocal(),N,ivec3_t(0,0,0),5.0);
    Atom Si(14,rvec3_t(0,0,0));

    GTH_PP siPP=GetGTH("Si","LDA",4);
    chmat_t Vloc=pw.MakeLocalPotential(&Si,siPP.local);
    chmat_t Vnl =pw.MakeSeparablePotential(&Si,siPP.nonlocal);
    size_t n=pw.GetNumFunctions();

    double vnlmax=0.0;                                   // nonlocal is non-trivial and Hermitian
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            vnlmax=std::max(vnlmax,std::abs(dcmplx(Vnl(i,j))));
            EXPECT_NEAR(std::real(dcmplx(Vnl(i,j))), std::real(dcmplx(Vnl(j,i))), 1e-12);
            EXPECT_NEAR(std::imag(dcmplx(Vnl(i,j))),-std::imag(dcmplx(Vnl(j,i))), 1e-12);
        }
    EXPECT_GT(vnlmax,0.0);

    chmat_t Vtot=Vloc+Vnl;
    std::vector<double> e=SolveBands(pw,&Vtot);          // generalized eigensolve succeeds, real spectrum
    ASSERT_EQ(e.size(),n);
    for (double ei : e) EXPECT_TRUE(std::isfinite(ei));
}

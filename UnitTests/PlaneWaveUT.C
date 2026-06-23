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
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.LASolver;
import qchem.Types;
import qchem.Blaze;
import qchem.Math;           // Pi

using BasisSet::Lattice_3D::PlaneWave_IBS;
using BasisSet::Lattice_3D::LocalPotential;
using BasisSet::Lattice_3D::BareCoulomb;
using BasisSet::Lattice_3D::GaussianSmearedNucleus;

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

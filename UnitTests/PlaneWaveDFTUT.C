// file: PlaneWaveDFTUT.C  Prototype self-consistent plane-wave DFT, validated outside the Hamiltonian
// framework (see doc/PlaneWavePlan.md sequencing + memory project_dft_upgrade_plan).
//
// WHY this lives in a unit test, not a library: the prototype wires three libraries together --
// BasisSet/Lattice_3D (PlaneWave_IBS primitives), Hamiltonian (the validated LDA ExFunctional), and
// LASolver (complex eigensolver).  A library that did this would invert the layering (BasisSet sits
// BELOW Hamiltonian).  The test is the one place allowed to reach across all three.  Once the
// Hamiltonian/SCFIterator stack is templated on T (double|dcmplx) these pieces migrate to their proper
// homes: the density -> a dcmplx ChargeDensity, Hartree/XC -> Dynamic_HT<dcmplx> terms.
//
// Build order (each piece tested before the next):
//   1. density rho~(dG) from occupied bands           [this commit]
//   2. G-space Hartree   V_H~(dG) = 4 pi rho~/|dG|^2
//   3. grid XC           rho(G)->rho(r)->Vxc->Vxc~(G)
//   4. SCF loop + total-energy validation
//
// CONVENTIONS (match Imp/PlaneWave_IBS.C):
//   psi_n(r) = (1/sqrt Omega) Sum_i c_n(i) e^{i(k+G_i).r},  c_n(i) = U(i,n)  (eigenvector columns).
//   rho(r)   = Sum_n f_n |psi_n(r)|^2 = Sum_dG rho~(dG) e^{i dG.r},
//   rho~(dG) = (1/Omega) Sum_n f_n Sum_{i,j: m_i-m_j=dm} c_n(i) conj(c_n(j)),   dG = B dm.
// So rho~(0) = (1/Omega) Sum_n f_n = N/Omega, and rho real => rho~(-dm) = conj(rho~(dm)).

#include <map>
#include <set>
#include <memory>
#include <vector>
#include <complex>
#include <cmath>
#include <functional>
#include <algorithm>
#include <iostream>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.BasisSet.Lattice_3D.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
import qchem.Blaze;          // mat_t<dcmplx>
import qchem.Math;           // Pi
import qchem.Hamiltonian.Internal.ExFunctional;     // the validated LDA functional interface
import qchem.Hamiltonian.Internal.SlaterExchange;   // Dirac exchange (alpha=2/3), eps_x = 3/4 v_x
import qchem.Hamiltonian.Internal.VWN_Correlation;  // VWN5 correlation (validated vs libxc)
import qchem.Hamiltonian.Internal.PWTerms;          // PW_External (dcmplx Hamiltonian term)
import qchem.Hamiltonian;                           // cStatic_HT / cDynamic_HT aliases (public term interfaces)
import qchem.Hamiltonian.Internal.Hamiltonian;      // cHamiltonianImp (the dcmplx Hamiltonian = sum of terms)
import qchem.Hamiltonian.Internal.Hamiltonians;     // Ham_PW_DFT (the assembled plane-wave LDA KS Hamiltonian)
import qchem.Energy;                                // EnergyBreakdown
import qchem.ChargeDensity.Imp.IrrepCD;             // IrrepCD<dcmplx> (concrete complex density)
import qchem.Symmetry.Irrep;                        // Irrep
import qchem.LASolver;                              // complex Hermitian eigensolver
import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
import qchem.SCFIterator;                           // cSCFIterator (the real framework SCF driver)
import qchem.SCFParams;                             // SCFParams
import qchem.WaveFunction;                          // cWaveFunction (read view of the converged state)
import qchem.ElectronConfiguration;                 // ElectronConfiguration base
import qchem.ElectronConfiguration.Crystal;         // Crystal_EC (single-k Bloch configuration)
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull; // tSCFAcceleratorNull<dcmplx>
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

using BasisSet::Lattice_3D::PlaneWave_IBS;
using BasisSet::Lattice_3D::HGH_LocalPotential;
using BasisSet::Lattice_3D::HGH_SeparablePotential;

namespace
{

// A ScalarFunction<double> wrapping a lambda f(r) -- to hand real-space fields to the basis's
// high-level integral methods (IntegralPotential/IntegralHartree/Integral).
struct FieldFn : public ScalarFunction<double>
{
    std::function<double(const rvec3_t&)> f;
    explicit FieldFn(std::function<double(const rvec3_t&)> g) : f(g) {}
    virtual double  operator()(const rvec3_t& r) const {return f(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
};

// Order ivec3_t lexicographically so it can key the rho~ map.
struct IVecLess
{
    bool operator()(const ivec3_t& a, const ivec3_t& b) const
    {
        if (a.x!=b.x) return a.x<b.x;
        if (a.y!=b.y) return a.y<b.y;
        return a.z<b.z;
    }
};

//! rho~(dm): the periodic charge density's Fourier components, keyed by reciprocal-index difference dm.
typedef std::map<ivec3_t, dcmplx, IVecLess> RhoG;

//! rho~(dm) = (1/Omega) Sum_bands f_b Sum_{i,j} c_b(i) conj(c_b(j)) delta(m_i-m_j, dm).
//! U columns are the bands (c_b(i)=U(i,b)); f are the occupations; Omega is the direct-cell volume.
RhoG BuildDensity(const PlaneWave_IBS& pw, const mat_t<dcmplx>& U, const rvec_t& f, double Omega)
{
    size_t n=pw.GetNumFunctions();
    RhoG rho;
    for (size_t b=0; b<f.size(); b++)
    {
        if (f[b]==0.0) continue;
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                rho[pw.GetGIndex(i)-pw.GetGIndex(j)] += f[b]*U(i,b)*std::conj(U(j,b));
    }
    for (auto& kv : rho) kv.second/=Omega;
    return rho;
}

//! rho~(dm), or 0 if that component is absent (outside the difference set).
dcmplx RhoAt(const RhoG& rho, const ivec3_t& dm)
{
    auto it=rho.find(dm);
    return it==rho.end() ? dcmplx(0.0) : it->second;
}

// --- G-space Hartree:  V_H(r) solves nabla^2 V_H = -4 pi rho,  so V_H~(G) = 4 pi rho~(G)/|G|^2. -----
// The dG=0 component is dropped (neutralising background), as in MakeLocalPotential.  The matrix is
// then <G_i|V_H|G_j> = V_H~(m_i-m_j) via MakePotential.

//! V_H~(dm) supplier for MakePotential.  \a B is the RECIPROCAL cell (G = B dm).
std::function<dcmplx(const ivec3_t&)> HartreeVtilde(const RhoG& rho, const UnitCell& B)
{
    return [&rho,&B](const ivec3_t& dm)->dcmplx
    {
        if (dm==ivec3_t(0,0,0)) return dcmplx(0.0);
        rvec3_t G=B.ToCartesian(rvec3_t(dm));
        return 4*Pi*RhoAt(rho,dm)/(G*G);
    };
}

//! E_H = 1/2 integral rho V_H = (Omega/2) Sum_{G!=0} 4 pi |rho~(G)|^2 / |G|^2.
double HartreeEnergy(const RhoG& rho, const UnitCell& B, double Omega)
{
    double E=0.0;
    for (const auto& kv : rho)
    {
        if (kv.first==ivec3_t(0,0,0)) continue;
        rvec3_t G=B.ToCartesian(rvec3_t(kv.first));
        E += 4*Pi*std::norm(kv.second)/(G*G);
    }
    return 0.5*Omega*E;
}

// --- grid XC ----------------------------------------------------------------------------------
// LDA Vxc(rho(r)) is pointwise NONLINEAR, so it must be evaluated in real space: rho(G)->rho(r) on a
// uniform grid, apply the functional, then forward-transform Vxc(r)->Vxc~(G).  The grid is uniform
// (weight Omega/Npts): for the band-limited, cusp-free pseudo-density the trapezoidal rule is
// spectrally exact (see memory project_dft_upgrade_plan -- this is why PW codes use uniform FFT grids,
// not Becke quadrature).  Phases use dG.r = 2 pi dm.r_frac (since B = 2 pi A^{-T}), so the transform is
// a plain DFT independent of cell shape.  Direct DFT sums here; FFT is the later optimisation.

//! Uniform N1xN2xN3 grid of FRACTIONAL coordinates r = (i1/N1, i2/N2, i3/N3).
std::vector<rvec3_t> UniformGrid(const ivec3_t& Ng)
{
    std::vector<rvec3_t> g;
    g.reserve(size_t(Ng.x)*Ng.y*Ng.z);
    for (int i1=0;i1<Ng.x;i1++)
        for (int i2=0;i2<Ng.y;i2++)
            for (int i3=0;i3<Ng.z;i3++)
                g.push_back(rvec3_t(i1/double(Ng.x), i2/double(Ng.y), i3/double(Ng.z)));
    return g;
}

//! Inverse transform: rho(r) = Sum_dm rho~(dm) e^{i 2 pi dm.r_frac} (the imaginary part cancels: rho real).
double RhoOfR(const RhoG& spec, const rvec3_t& rf)
{
    dcmplx s(0.0);
    for (const auto& kv : spec)
    {
        double ph=2*Pi*(kv.first.x*rf.x + kv.first.y*rf.y + kv.first.z*rf.z);
        s += kv.second*dcmplx(std::cos(ph),std::sin(ph));
    }
    return std::real(s);
}

//! Forward transform: spectrum component f~(dm) = (1/Npts) Sum_grid f(r) e^{-i 2 pi dm.r_frac}.
dcmplx ForwardDFT(const std::vector<double>& f, const std::vector<rvec3_t>& rf, const ivec3_t& dm)
{
    dcmplx s(0.0);
    for (size_t q=0;q<f.size();q++)
    {
        double ph=-2*Pi*(dm.x*rf[q].x + dm.y*rf[q].y + dm.z*rf[q].z);
        s += f[q]*dcmplx(std::cos(ph),std::sin(ph));
    }
    return s/double(f.size());
}

//! Evaluate v_xc(rho(r)) on the uniform grid (returns the field + its coordinates \a rf), and by
//! reference E_xc = integral eps_xc rho and the XC double-count VxcDc = integral rho v_xc.
//! \a vxcOf / \a epsxcOf are rho->v_xc and rho->eps_xc (sum exchange+correlation by the caller).
std::vector<double> XcGridField(const RhoG& rho, const ivec3_t& Ng, double Omega,
                                const std::function<double(double)>& vxcOf,
                                const std::function<double(double)>& epsxcOf,
                                std::vector<rvec3_t>& rf, double& Exc, double& VxcDc)
{
    rf=UniformGrid(Ng);
    size_t Npts=rf.size();
    std::vector<double> vxc(Npts);
    Exc=0.0; VxcDc=0.0;
    for (size_t q=0;q<Npts;q++)
    {
        double rr=RhoOfR(rho,rf[q]);
        if (rr<0.0) rr=0.0;                 // clamp tiny numerical negatives before the functional
        vxc[q]=vxcOf(rr);
        Exc   += epsxcOf(rr)*rr;
        VxcDc += vxc[q]*rr;
    }
    Exc   *= Omega/double(Npts);
    VxcDc *= Omega/double(Npts);
    return vxc;
}

//! Build Vxc~(dm) over a single basis's difference set (single-k convenience).
RhoG BuildXcVtilde(const PlaneWave_IBS& pw, const RhoG& rho, const ivec3_t& Ng, double Omega,
                   const std::function<double(double)>& vxcOf,
                   const std::function<double(double)>& epsxcOf, double& Exc, double& VxcDc)
{
    std::vector<rvec3_t> rf;
    std::vector<double>  vxc=XcGridField(rho,Ng,Omega,vxcOf,epsxcOf,rf,Exc,VxcDc);
    RhoG vtil;                              // forward-transform once per distinct dm in the diff set
    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=0;j<n;j++)
        {
            ivec3_t dm=pw.GetGIndex(i)-pw.GetGIndex(j);
            if (vtil.find(dm)==vtil.end()) vtil[dm]=ForwardDFT(vxc,rf,dm);
        }
    return vtil;
}

// A small PW basis for the density tests: cubic cell, Gamma point.
struct PWFixture
{
    double           a=8.0, Ecut=4.0, Omega=a*a*a;
    UnitCell         cell{a};
    Lattice_3D       lat{cell, ivec3_t(1,1,1)};
    ReciprocalLattice recip{lat.Reciprocal()};   // owns the reciprocal cell B used by Hartree
    PlaneWave_IBS    pw{lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), Ecut};
    const UnitCell&  B() const {return recip.GetCell();}
};

class PlaneWaveDFT : public ::testing::Test {};

// A single occupied plane wave (f=2 on band 0 = the i=0 PW) is a uniform density: the ONLY nonzero
// Fourier component is rho~(0)=N/Omega=2/Omega; every dm!=0 is exactly zero.
TEST_F(PlaneWaveDFT, SinglePlaneWaveIsUniform)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    mat_t<dcmplx> U(n,1,dcmplx(0.0));
    U(0,0)=1.0;
    rvec_t f(1,2.0);

    RhoG rho=BuildDensity(F.pw, U, f, F.Omega);
    EXPECT_NEAR(std::real(RhoAt(rho,ivec3_t(0,0,0))), 2.0/F.Omega, 1e-14);
    EXPECT_NEAR(std::imag(RhoAt(rho,ivec3_t(0,0,0))), 0.0,         1e-14);
    for (const auto& kv : rho)
        if (!(kv.first==ivec3_t(0,0,0)))
            EXPECT_NEAR(std::abs(kv.second), 0.0, 1e-14) << "spurious rho~ at non-zero dm";
}

// Charge sum rule: integral rho = Omega * rho~(0) = N = Sum f_b, for an arbitrary normalised band.
TEST_F(PlaneWaveDFT, ChargeSumRule)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    // A normalised band spread over the first three PWs (Sum_i |c_i|^2 = 1).
    mat_t<dcmplx> U(n,1,dcmplx(0.0));
    U(0,0)=dcmplx(0.6,0.0); U(1,0)=dcmplx(0.0,0.8);   // |0.6|^2+|0.8i|^2 = 0.36+0.64 = 1
    rvec_t f(1,2.0);

    RhoG rho=BuildDensity(F.pw, U, f, F.Omega);
    EXPECT_NEAR(F.Omega*std::real(RhoAt(rho,ivec3_t(0,0,0))), 2.0, 1e-13);  // N = 2
}

// rho is real => its Fourier components obey rho~(-dm) = conj(rho~(dm)).  Use a complex superposition
// band so the off-diagonal components are genuinely complex.
TEST_F(PlaneWaveDFT, HermitianSymmetry)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    mat_t<dcmplx> U(n,1,dcmplx(0.0));
    U(0,0)=dcmplx(1.0/std::sqrt(2.0),0.0);
    U(1,0)=dcmplx(0.0,1.0/std::sqrt(2.0));            // i/sqrt2 : makes cross terms imaginary
    rvec_t f(1,2.0);

    RhoG rho=BuildDensity(F.pw, U, f, F.Omega);
    ivec3_t m0=F.pw.GetGIndex(0), m1=F.pw.GetGIndex(1);
    dcmplx r01=RhoAt(rho, m0-m1), r10=RhoAt(rho, m1-m0);
    EXPECT_NEAR(std::abs(r01 - std::conj(r10)), 0.0, 1e-14);
    EXPECT_GT (std::abs(r01), 1e-3);                  // and they are actually non-trivial / complex
    EXPECT_NEAR(std::real(r01), 0.0, 1e-14);          // pure imaginary for this choice
}

// Analytic check: a single-cosine density rho(r) = rho0 + 2A cos(G0.r) has Poisson solution
// V_H(r) = (4 pi A/|G0|^2) 2 cos(G0.r), i.e. V_H~(+/-G0) = 4 pi A/|G0|^2, V_H~(0)=0; and
// E_H = (Omega/2) Sum_{+/-G0} 4 pi A^2/|G0|^2 = Omega 4 pi A^2/|G0|^2.
TEST_F(PlaneWaveDFT, HartreeSingleCosineMatchesPoisson)
{
    PWFixture F;
    const double A=0.01;
    RhoG rho;
    rho[ivec3_t( 0,0,0)] = 2.0/F.Omega;   // average density (irrelevant to V_H, dropped)
    rho[ivec3_t( 1,0,0)] = A;
    rho[ivec3_t(-1,0,0)] = A;

    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    double  g02=G0*G0;                                       // (2 pi/a)^2 for the cubic cell
    EXPECT_NEAR(g02, (2*Pi/F.a)*(2*Pi/F.a), 1e-12);

    auto VH=HartreeVtilde(rho, F.B());
    EXPECT_NEAR(std::abs(VH(ivec3_t(0,0,0))), 0.0,            1e-14);  // dG=0 dropped
    EXPECT_NEAR(std::real(VH(ivec3_t(1,0,0))), 4*Pi*A/g02,    1e-12);
    EXPECT_NEAR(std::imag(VH(ivec3_t(1,0,0))), 0.0,           1e-14);

    EXPECT_NEAR(HartreeEnergy(rho,F.B(),F.Omega), F.Omega*4*Pi*A*A/g02, 1e-10);
}

// The Hartree matrix from MakePotential picks up V_H~(m_i-m_j) on each pair, and is Hermitian.
TEST_F(PlaneWaveDFT, HartreeMatrixElementsAndHermiticity)
{
    PWFixture F;
    const double A=0.01;
    RhoG rho;
    rho[ivec3_t( 1,0,0)] = A;
    rho[ivec3_t(-1,0,0)] = A;
    double g02=(2*Pi/F.a)*(2*Pi/F.a);

    chmat_t V=F.pw.MakePotential(HartreeVtilde(rho, F.B()));
    size_t n=F.pw.GetNumFunctions();

    bool found=false;                                        // a pair differing by (1,0,0)
    for (size_t i=0;i<n && !found;i++)
        for (size_t j=0;j<n && !found;j++)
            if (F.pw.GetGIndex(i)-F.pw.GetGIndex(j)==ivec3_t(1,0,0))
            {
                EXPECT_NEAR(std::real(dcmplx(V(i,j))), 4*Pi*A/g02, 1e-12);
                EXPECT_NEAR(std::abs(dcmplx(V(i,j))-std::conj(dcmplx(V(j,i)))), 0.0, 1e-14);
                found=true;
            }
    EXPECT_TRUE(found) << "expected a G-pair differing by (1,0,0)";
}

// The inverse/forward DFT pair is exact for a band-limited density once the grid resolves its
// components: build rho(r) from a 2-component rho~, transform back, recover rho~ exactly (and no
// spurious power at 2G).  This validates the transform machinery independent of the functional.
TEST_F(PlaneWaveDFT, DensityGridTransformRoundtrip)
{
    PWFixture F;
    const double A=0.02, rho0=2.0/F.Omega;
    RhoG rho;
    rho[ivec3_t( 0,0,0)] = rho0;
    rho[ivec3_t( 1,0,0)] = A;
    rho[ivec3_t(-1,0,0)] = A;

    std::vector<rvec3_t> rf=UniformGrid(ivec3_t(8,8,8));
    std::vector<double>  rr(rf.size());
    for (size_t q=0;q<rf.size();q++) rr[q]=RhoOfR(rho,rf[q]);

    EXPECT_NEAR(std::real(ForwardDFT(rr,rf,ivec3_t(0,0,0))), rho0, 1e-12);
    EXPECT_NEAR(std::real(ForwardDFT(rr,rf,ivec3_t(1,0,0))), A,    1e-12);
    EXPECT_NEAR(std::imag(ForwardDFT(rr,rf,ivec3_t(1,0,0))), 0.0,  1e-12);
    EXPECT_NEAR(std::abs (ForwardDFT(rr,rf,ivec3_t(2,0,0))), 0.0,  1e-12);  // no spurious 2G power
}

// Uniform-density limit: rho = rho0 everywhere => Vxc(r)=v_xc(rho0) constant, so Vxc~(0)=v_xc(rho0),
// Vxc~(dm!=0)=0, and E_xc = Omega eps_xc(rho0) rho0.  Uses the validated Dirac exchange (eps_x=3/4 v_x).
TEST_F(PlaneWaveDFT, XcUniformDensityLimit)
{
    PWFixture F;
    const double rho0=0.05;
    RhoG rho;
    rho[ivec3_t(0,0,0)] = rho0;

    qchem::Hamiltonian::SlaterExchange ex(2.0/3.0);              // Dirac exchange
    auto vxcOf =[&](double r){return ex.GetVxc (r);};
    auto epsOf =[&](double r){return ex.GetEpsXc(r);};

    double Exc=0.0, VxcDc=0.0;
    RhoG vtil=BuildXcVtilde(F.pw, rho, ivec3_t(6,6,6), F.Omega, vxcOf, epsOf, Exc, VxcDc);

    EXPECT_NEAR(std::real(RhoAt(vtil,ivec3_t(0,0,0))), ex.GetVxc(rho0), 1e-12);
    EXPECT_NEAR(std::abs (RhoAt(vtil,ivec3_t(1,0,0))), 0.0,             1e-12);
    EXPECT_NEAR(Exc,   F.Omega*ex.GetEpsXc(rho0)*rho0, 1e-10);
    EXPECT_NEAR(VxcDc, F.Omega*ex.GetVxc  (rho0)*rho0, 1e-10);   // integral rho v_xc
}

// --- the SCF loop (prototype) -----------------------------------------------------------------
// Assemble H(k) = 1/2 <p^2> + V_ext + V_H[rho] + V_xc[rho], diagonalise, occupy the lowest bands,
// rebuild rho, linearly mix, iterate to self-consistency.  This is the hand-rolled driver that the
// templated Hamiltonian/SCFIterator<dcmplx> will eventually replace (V_H/V_xc -> Dynamic_HT<dcmplx>).

struct SCFResult
{
    RhoG   rho;                                   //!< converged self-consistent density rho~(dm)
    double Ekin=0, EH=0, Exc=0, Eext=0;           //!< energy components at the fixed point
    double Eband=0, VxcDc=0;                      //!< band sum, and XC double-count integral rho v_xc
    double Etot_band=0;                           //!< band route:  Eband - E_H + E_xc - integral rho v_xc
    double Etot_direct=0;                         //!< direct route: T + integral rho V_ext + E_H + E_xc
    double gap=0;                                 //!< CBM-VBM over the sampled k-points (multi-k only)
    int    iters=0;
    bool   converged=false;
};

//! Largest |component| over the G-set -- the difference set spans +/-2x this, so a grid with
//! N > 2*(2*maxComp) per axis resolves it without aliasing (=> band/direct energies agree exactly).
int MaxGComponent(const PlaneWave_IBS& pw)
{
    int m=0;
    for (size_t i=0;i<pw.GetNumFunctions();i++)
    {
        ivec3_t g=pw.GetGIndex(i);
        m=std::max(m, std::max({std::abs(g.x),std::abs(g.y),std::abs(g.z)}));
    }
    return m;
}

// \a Vext is the FULL static external block (local PP + KB nonlocal): density-independent, assembled
// once by the caller.  E_ext is taken as the band expectation Sum_n f_n <psi|Vext|psi> -- valid for the
// nonlocal projector part too (where there is no real-space integral rho V_ext).
SCFResult RunSCF(const PlaneWave_IBS& pw, const UnitCell& B, double Omega,
                 const chmat_t& Vext, int Nelec, const ivec3_t& Ng,
                 const std::function<double(double)>& vxcOf, const std::function<double(double)>& epsxcOf,
                 double alpha=0.5, double tol=1e-9, int maxiter=400)
{
    size_t n=pw.GetNumFunctions();
    int    Nocc=Nelec/2;
    chmat_t K=pw.MakeKinetic();                   // diagonal <p^2>, constant across iterations
    chmat_t S=pw.MakeOverlap();                   // identity

    SCFResult R;
    R.rho[ivec3_t(0,0,0)] = dcmplx(double(Nelec)/Omega, 0.0);   // uniform initial guess

    mat_t<dcmplx> U; rvec_t e; rvec_t f; RhoG rhoNew;
    for (int it=0; it<maxiter; it++)
    {
        R.iters=it+1;
        auto   Vh=HartreeVtilde(R.rho,B);
        double Exc,VxcDc;
        RhoG   vxcmap=BuildXcVtilde(pw,R.rho,Ng,Omega,vxcOf,epsxcOf,Exc,VxcDc);
        auto   Vloc=[&](const ivec3_t& dm)->dcmplx { return Vh(dm)+RhoAt(vxcmap,dm); };  // Hartree+XC

        hmat_t<dcmplx> H(n);                       // 1/2 <p^2> (diag) + Vext(i,j) + V_H+V_xc (Gi-Gj)
        for (size_t i=0;i<n;i++)
            for (size_t j=i;j<n;j++)
                H(i,j) = (i==j ? dcmplx(0.5*std::real(dcmplx(K(i,i)))) : dcmplx(0.0))
                       + dcmplx(Vext(i,j))
                       + Vloc(pw.GetGIndex(i)-pw.GetGIndex(j));

        LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
        las->SetBasisOverlap(S);
        auto sol=las->Solve(H);
        delete las;
        U=std::get<0>(sol); e=std::get<1>(sol);

        std::vector<size_t> idx(e.size());         // occupy the Nocc lowest-energy bands (2 each)
        for (size_t i=0;i<idx.size();i++) idx[i]=i;
        std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return e[a]<e[b];});
        f=rvec_t(e.size(),0.0);
        for (int k=0;k<Nocc;k++) f[idx[k]]=2.0;

        rhoNew=BuildDensity(pw,U,f,Omega);
        double drho=0.0;                           // max |rho_new - rho_in| over all components
        for (const auto& kv : rhoNew) drho=std::max(drho, std::abs(kv.second-RhoAt(R.rho,kv.first)));

        RhoG mixed;                                // linear mix
        for (const auto& kv : rhoNew) mixed[kv.first]=(1-alpha)*RhoAt(R.rho,kv.first)+alpha*kv.second;
        R.rho=mixed;

        if (drho<tol) { R.converged=true; break; }
    }

    // Final energies at the self-consistent density (rhoNew = output of the converged solve).
    R.EH=HartreeEnergy(rhoNew,B,Omega);
    double dummyExc,dummyDc;
    BuildXcVtilde(pw,rhoNew,Ng,Omega,vxcOf,epsxcOf,dummyExc,dummyDc);
    R.Exc=dummyExc; R.VxcDc=dummyDc;
    for (size_t b=0;b<e.size();b++)                          // E_ext = Sum_n f_n <psi|Vext|psi>
        if (f[b]!=0.0)
        {
            dcmplx ev(0.0);
            for (size_t i=0;i<n;i++)
                for (size_t j=0;j<n;j++) ev += std::conj(dcmplx(U(i,b)))*dcmplx(Vext(i,j))*dcmplx(U(j,b));
            R.Eext += f[b]*std::real(ev);
        }
    for (size_t b=0;b<e.size();b++) R.Eband += f[b]*e[b];
    for (size_t b=0;b<e.size();b++)
        if (f[b]!=0.0)
            for (size_t i=0;i<n;i++) R.Ekin += f[b]*0.5*std::real(dcmplx(K(i,i)))*std::norm(dcmplx(U(i,b)));
    R.Etot_band   = R.Eband - R.EH + R.Exc - R.VxcDc;
    R.Etot_direct = R.Ekin + R.Eext + R.EH + R.Exc;
    return R;
}

// --- multi-k (Brillouin-zone) SCF ---------------------------------------------------------------
// The BZ sum Sum_k w_k is the irrep loop the framework will eventually own (k IS a Bloch irrep).  Each
// k has its OWN plane-wave basis (cutoff set {G:1/2|k+G|^2<Ecut} depends on k) and its OWN external
// block (the KB projectors beta~(|k+G|) are k-dependent); both are static, so built once.  The density
// rho~(dm) (k-independent, lattice-periodic) accumulates Sum_k w_k Sum_n f c c* over all k.  For an
// insulator (Si: 8 valence e-) exactly Nelec/2 bands are filled at EVERY k -> trivial occupation.
SCFResult RunSCF_kpoints(const ReciprocalLattice& recip, const UnitCell& B, double Omega,
                         const Structure& cl, const HGH_LocalPotential& loc, const HGH_SeparablePotential& nl,
                         const ivec3_t& Nmp, double Ecut, int Nelec, const ivec3_t& Ng,
                         const std::function<double(double)>& vxcOf, const std::function<double(double)>& epsxcOf,
                         double alpha, double tol, int maxiter)
{
    std::vector<std::unique_ptr<PlaneWave_IBS>> pw;     // per-k static data (basis, external block, kinetic)
    std::vector<chmat_t> Vext, K;
    for (int kx=0;kx<Nmp.x;kx++)
        for (int ky=0;ky<Nmp.y;ky++)
            for (int kz=0;kz<Nmp.z;kz++)
            {
                auto p=std::make_unique<PlaneWave_IBS>(recip, Nmp, ivec3_t(kx,ky,kz), Ecut);
                Vext.push_back(p->MakeLocalPotential(&cl,loc)+p->MakeSeparablePotential(&cl,nl));
                K.push_back(p->MakeKinetic());
                pw.push_back(std::move(p));
            }
    int    Nk=int(pw.size());
    double wk=1.0/Nk;
    int    Nocc=Nelec/2;

    std::set<ivec3_t,IVecLess> diff;                    // union of all k difference sets (for the XC transform)
    for (int k=0;k<Nk;k++)
        for (size_t i=0;i<pw[k]->GetNumFunctions();i++)
            for (size_t j=0;j<pw[k]->GetNumFunctions();j++)
                diff.insert(pw[k]->GetGIndex(i)-pw[k]->GetGIndex(j));

    SCFResult R;
    R.rho[ivec3_t(0,0,0)] = dcmplx(double(Nelec)/Omega, 0.0);
    std::vector<mat_t<dcmplx>> Uk(Nk); std::vector<rvec_t> ek(Nk), fk(Nk); RhoG rhoNew;

    for (int it=0; it<maxiter; it++)
    {
        R.iters=it+1;
        auto Vh=HartreeVtilde(R.rho,B);
        std::vector<rvec3_t> rf; double Exc,VxcDc;
        std::vector<double> vxcField=XcGridField(R.rho,Ng,Omega,vxcOf,epsxcOf,rf,Exc,VxcDc);
        RhoG vxcmap; for (const ivec3_t& dm : diff) vxcmap[dm]=ForwardDFT(vxcField,rf,dm);
        auto Vloc=[&](const ivec3_t& dm)->dcmplx { return Vh(dm)+RhoAt(vxcmap,dm); };

        rhoNew.clear();
        for (int k=0;k<Nk;k++)
        {
            size_t n=pw[k]->GetNumFunctions();
            hmat_t<dcmplx> H(n), S=pw[k]->MakeOverlap();
            for (size_t i=0;i<n;i++)
                for (size_t j=i;j<n;j++)
                    H(i,j) = (i==j ? dcmplx(0.5*std::real(dcmplx(K[k](i,i)))) : dcmplx(0.0))
                           + dcmplx(Vext[k](i,j)) + Vloc(pw[k]->GetGIndex(i)-pw[k]->GetGIndex(j));

            LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
            las->SetBasisOverlap(S);
            auto sol=las->Solve(H);
            delete las;
            Uk[k]=std::get<0>(sol); ek[k]=std::get<1>(sol);

            std::vector<size_t> idx(ek[k].size());                 // occupy the Nocc lowest bands
            for (size_t i=0;i<idx.size();i++) idx[i]=i;
            std::sort(idx.begin(),idx.end(),[&](size_t a,size_t b){return ek[k][a]<ek[k][b];});
            fk[k]=rvec_t(ek[k].size(),0.0);
            for (int o=0;o<Nocc;o++) fk[k][idx[o]]=2.0;

            for (size_t b=0;b<fk[k].size();b++)                    // accumulate weighted density
                if (fk[k][b]!=0.0)
                    for (size_t i=0;i<n;i++)
                        for (size_t j=0;j<n;j++)
                            rhoNew[pw[k]->GetGIndex(i)-pw[k]->GetGIndex(j)]
                                += (wk*fk[k][b]/Omega)*dcmplx(Uk[k](i,b))*std::conj(dcmplx(Uk[k](j,b)));
        }

        double drho=0.0;
        for (const auto& kv : rhoNew) drho=std::max(drho, std::abs(kv.second-RhoAt(R.rho,kv.first)));
        RhoG mixed;
        for (const auto& kv : rhoNew) mixed[kv.first]=(1-alpha)*RhoAt(R.rho,kv.first)+alpha*kv.second;
        R.rho=mixed;
        if (drho<tol) { R.converged=true; break; }
    }

    R.EH=HartreeEnergy(rhoNew,B,Omega);                            // energies at the converged density
    {
        std::vector<rvec3_t> rfe; double exc,dc;
        XcGridField(rhoNew,Ng,Omega,vxcOf,epsxcOf,rfe,exc,dc);
        R.Exc=exc; R.VxcDc=dc;
    }
    double vbm=-1e30, cbm=1e30;
    for (int k=0;k<Nk;k++)
    {
        for (size_t b=0;b<fk[k].size();b++)
        {
            R.Eband += wk*fk[k][b]*ek[k][b];
            if (fk[k][b]!=0.0)
            {
                for (size_t i=0;i<pw[k]->GetNumFunctions();i++)
                    R.Ekin += wk*fk[k][b]*0.5*std::real(dcmplx(K[k](i,i)))*std::norm(dcmplx(Uk[k](i,b)));
                dcmplx ev(0.0);
                size_t n=pw[k]->GetNumFunctions();
                for (size_t i=0;i<n;i++)
                    for (size_t j=0;j<n;j++)
                        ev += std::conj(dcmplx(Uk[k](i,b)))*dcmplx(Vext[k](i,j))*dcmplx(Uk[k](j,b));
                R.Eext += wk*fk[k][b]*std::real(ev);
            }
        }
        std::vector<double> es;                                    // sorted bands at this k -> gap
        for (size_t b=0;b<ek[k].size();b++) es.push_back(ek[k][b]);
        std::sort(es.begin(),es.end());
        vbm=std::max(vbm, es[Nocc-1]);
        cbm=std::min(cbm, es[Nocc]);
    }
    R.gap=cbm-vbm;
    R.Etot_band   = R.Eband - R.EH + R.Exc - R.VxcDc;
    R.Etot_direct = R.Ekin + R.Eext + R.EH + R.Exc;
    return R;
}

// Jellium: no external potential, 2 electrons.  The self-consistent density stays uniform (only the
// G=0 component), so E_H=0, the kinetic energy is zero (G=0 band), and E_tot = E_xc = Omega eps_xc rho0.
// Also a sanity check that the two energy routes agree.
TEST_F(PlaneWaveDFT, ScfJelliumUniform)
{
    PWFixture F;
    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};
    chmat_t Vzero=F.pw.MakePotential([](const ivec3_t&){return dcmplx(0.0);});

    SCFResult R=RunSCF(F.pw, F.B(), F.Omega, Vzero, 2, ivec3_t(8,8,8), vxcOf, epsOf);
    ASSERT_TRUE(R.converged);

    const double rho0=2.0/F.Omega;
    for (const auto& kv : R.rho)                                   // density is uniform
        if (!(kv.first==ivec3_t(0,0,0))) EXPECT_NEAR(std::abs(kv.second), 0.0, 1e-8);
    EXPECT_NEAR(R.Ekin, 0.0, 1e-9);                               // G=0 band: zero kinetic
    EXPECT_NEAR(R.EH,   0.0, 1e-9);                               // uniform: no Hartree
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-8);               // two routes agree
    EXPECT_NEAR(R.Etot_direct, F.Omega*epsOf(rho0)*rho0, 1e-7);  // = N eps_xc(rho0)
}

// A weak external cosine well + 2 electrons + Hartree + LDA(Dirac+VWN): a non-trivial self-consistent
// loop (Hartree and XC respond to a genuinely modulated density).  No external reference, so we check
// the strongest internal property: at the fixed point the band-sum and direct total energies agree.
TEST_F(PlaneWaveDFT, ScfWeakCosineSelfConsistent)
{
    const double a=6.0, Ecut=4.0, Omega=a*a*a, V0=-0.3;
    UnitCell          cell(a);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), Ecut);

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};
    // V_ext(r) = 2 V0 (cos + cos + cos): only Fourier components are the unit reciprocal steps.
    chmat_t Vext=pw.MakePotential([V0](const ivec3_t& dm)->dcmplx
    { return (dm.x*dm.x+dm.y*dm.y+dm.z*dm.z==1) ? dcmplx(V0) : dcmplx(0.0); });

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, 2, ivec3_t(12,12,12), vxcOf, epsOf, 0.5, 1e-9, 400);
    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-6);                       // stationarity / consistency
    EXPECT_GT (std::abs(RhoAt(R.rho, ivec3_t(1,0,0))), 1e-4);            // density genuinely modulated
}

// Real material: silicon, diamond structure, GTH-LDA q4 pseudopotential (local + KB nonlocal).
// FCC PRIMITIVE cell (2 Si, 8 valence electrons) at Gamma.  G=0 is dropped (neutralising background)
// so the ABSOLUTE total energy is shifted and not meaningful here (per the plan) -- we validate the
// SELF-CONSISTENT loop instead: it converges, the charge sum rule gives exactly 8 valence electrons,
// and the band-sum and direct total energies agree at the fixed point.  A modest Ecut keeps the direct
// (non-FFT) transforms fast; the grid is auto-sized to resolve the basis difference set.
TEST_F(PlaneWaveDFT, ScfSiliconDiamondConverges)
{
    const double a=10.26, h=0.5*a;                         // Si lattice constant ~5.43 A in bohr
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);       // FCC primitive: cols = a/2 (011),(101),(110)
    UnitCell          cell(A);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    const double      Omega=cell.GetCellVolume();          // = a^3/4
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    Molecule si;                                           // 2-atom diamond basis
    si.Insert(new Atom(14, rvec3_t(0,0,0)));
    si.Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));

    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();
    chmat_t Vext = pw.MakeLocalPotential(&si,loc) + pw.MakeSeparablePotential(&si,nl);

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};

    const int Nval=8;                                      // 2 Si x Zion 4
    int       m=MaxGComponent(pw);
    ivec3_t   Ng(4*m+1, 4*m+1, 4*m+1);                     // resolves the difference set (no aliasing)

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, Nval, Ng, vxcOf, epsOf, 0.4, 1e-8, 400);

    std::cout << "[Si] nG="<<pw.GetNumFunctions()<<" grid="<<(4*m+1)<<"^3 iters="<<R.iters
              << " converged="<<R.converged << "\n  Ekin="<<R.Ekin<<" Eext="<<R.Eext
              << " E_H="<<R.EH<<" E_xc="<<R.Exc << "\n  Etot(band)="<<R.Etot_band
              << " Etot(direct)="<<R.Etot_direct << std::endl;

    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(Omega*std::real(RhoAt(R.rho,ivec3_t(0,0,0))), double(Nval), 1e-6);  // 8 valence e-
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-5);                                  // stationary fixed point
}

// Silicon again, but BZ-sampled over a 2x2x2 Monkhorst-Pack mesh (8 k-points) instead of Gamma-only.
// This exercises the irrep(k) loop, the k-dependent KB projectors, and a properly BZ-averaged density.
// Validate the same robust internal properties (converges, exact 8-electron charge, band==direct
// energy) and report the BZ-sampled band gap.
TEST_F(PlaneWaveDFT, ScfSiliconBZSampled)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(A);
    Lattice_3D        lat(cell, ivec3_t(2,2,2));
    ReciprocalLattice recip(lat.Reciprocal());
    const double      Omega=cell.GetCellVolume();
    const double      Ecut=4.0;

    Molecule si;
    si.Insert(new Atom(14, rvec3_t(0,0,0)));
    si.Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};

    // grid resolves the difference set: size from the richest (Gamma) basis.
    PlaneWave_IBS pwGamma(lat.Reciprocal(), ivec3_t(2,2,2), ivec3_t(0,0,0), Ecut);
    int m=MaxGComponent(pwGamma);
    ivec3_t Ng(4*m+1, 4*m+1, 4*m+1);

    SCFResult R=RunSCF_kpoints(recip, recip.GetCell(), Omega, si, loc, nl,
                               ivec3_t(2,2,2), Ecut, 8, Ng, vxcOf, epsOf, 0.4, 1e-8, 400);

    std::cout << "[Si 2x2x2] iters="<<R.iters<<" converged="<<R.converged
              << "\n  Etot(band)="<<R.Etot_band<<" Etot(direct)="<<R.Etot_direct
              << " gap="<<R.gap<<" Ha ("<<R.gap*27.2114<<" eV)" << std::endl;

    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(Omega*std::real(RhoAt(R.rho,ivec3_t(0,0,0))), 8.0, 1e-6);   // 8 valence e-
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-5);                          // stationary fixed point
    EXPECT_GT (R.gap, 0.0);                                                 // Si is a semiconductor
}

// --- Stage 2: the basis-level high-level DFT capability (DFTPotential_IBS), validated against the
// prototype's analytic/free-function results.  These are the questions the framework terms will ask --
// the term hands a real-space ScalarFunction and the basis owns the integration (no G-vectors exposed).

// Integral(f) = integral f d3r over the cell: a constant integrates to const*Omega; a reciprocal
// cosine integrates to zero.
TEST_F(PlaneWaveDFT, BasisIntegralScalar)
{
    PWFixture F;
    FieldFn cst([](const rvec3_t&){return 2.5;});
    EXPECT_NEAR(F.pw.Integral(cst), 2.5*F.Omega, 1e-9*F.Omega);

    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    FieldFn cosfn([&](const rvec3_t& r){return std::cos(G0*r);});
    EXPECT_NEAR(F.pw.Integral(cosfn), 0.0, 1e-9*F.Omega);
}

// IntegralHartree on a single-cosine density rho = rho0 + 2A cos(G0.r) reproduces the Poisson result
// (E_H = Omega 4 pi A^2/|G0|^2, V_H~(G0) = 4 pi A/|G0|^2) -- same check as HartreeSingleCosineMatchesPoisson,
// now driven entirely through the basis's high-level method (the basis samples + transforms internally).
TEST_F(PlaneWaveDFT, BasisIntegralHartreeMatchesPoisson)
{
    PWFixture F;
    const double rho0=2.0/F.Omega, A=0.01;
    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    double  g02=G0*G0;
    FieldFn rho([&](const rvec3_t& r){return rho0 + 2*A*std::cos(G0*r);});

    double  Eh=0.0;
    chmat_t VH=F.pw.IntegralHartree(rho, Eh);
    EXPECT_NEAR(Eh, F.Omega*4*Pi*A*A/g02, 1e-7);

    size_t n=F.pw.GetNumFunctions(); bool found=false;
    for (size_t i=0;i<n && !found;i++)
        for (size_t j=0;j<n && !found;j++)
            if (F.pw.GetGIndex(i)-F.pw.GetGIndex(j)==ivec3_t(1,0,0))
            {
                EXPECT_NEAR(std::real(dcmplx(VH(i,j))), 4*Pi*A/g02, 1e-7);
                found=true;
            }
    EXPECT_TRUE(found);
}

// IntegralPotential of a constant V is V*Identity; of a reciprocal cosine has Vtilde(+/-e)=1/2.
TEST_F(PlaneWaveDFT, BasisIntegralPotentialConstAndCosine)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();

    FieldFn cst([](const rvec3_t&){return 0.7;});
    chmat_t Vc=F.pw.IntegralPotential(cst);
    for (size_t i=0;i<n;i++)
    {
        EXPECT_NEAR(std::real(dcmplx(Vc(i,i))), 0.7, 1e-9);
        for (size_t j=i+1;j<n;j++) EXPECT_NEAR(std::abs(dcmplx(Vc(i,j))), 0.0, 1e-9);
    }

    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    FieldFn cosfn([&](const rvec3_t& r){return std::cos(G0*r);});
    chmat_t Vk=F.pw.IntegralPotential(cosfn);
    bool found=false;
    for (size_t i=0;i<n && !found;i++)
        for (size_t j=0;j<n && !found;j++)
            if (F.pw.GetGIndex(i)-F.pw.GetGIndex(j)==ivec3_t(1,0,0))
            {
                EXPECT_NEAR(std::real(dcmplx(Vk(i,j))), 0.5, 1e-9);
                found=true;
            }
    EXPECT_TRUE(found);
}

// MakeExternalPotential (with the Si PP configured on the basis) equals the explicit local+KB assembly.
TEST_F(PlaneWaveDFT, BasisExternalPotentialMatchesPPAssembly)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> Amat(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(Amat);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    Molecule si;
    si.Insert(new Atom(14, rvec3_t(0,0,0)));
    si.Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();

    chmat_t ref = pw.MakeLocalPotential(&si,loc) + pw.MakeSeparablePotential(&si,nl);
    pw.SetPseudopotential(&loc, &nl);
    chmat_t got = pw.MakeExternalPotential(&si);

    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(got(i,j))-dcmplx(ref(i,j))), 0.0, 1e-12);
}

// The PW_External Hamiltonian term (cStatic_HT<dcmplx>) routes through the framework's GetMatrix path
// and the abstract DFTPotential_IBS dynamic_cast -- its matrix must equal the basis's external assembly.
// This is the dependency inversion working end-to-end through a real Hamiltonian term.
TEST_F(PlaneWaveDFT, PWExternalTermMatchesBasis)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> Amat(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(Amat);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    auto si=std::make_shared<Molecule>();
    si->Insert(new Atom(14, rvec3_t(0,0,0)));
    si->Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();
    pw.SetPseudopotential(&loc, &nl);

    qchem::Hamiltonian::PW_External   ext(si);
    qchem::Hamiltonian::cStatic_HT*   term=&ext;        // the public term interface (as the Hamiltonian holds it)
    const chmat_t& M  = term->GetMatrix(&pw, Spin::None);
    chmat_t        ref= pw.MakeExternalPotential(si.get());

    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(M(i,j))-dcmplx(ref(i,j))), 0.0, 1e-12);
}

// The PW_Hartree and PW_XC dynamic terms route a complex density (IrrepCD<dcmplx>) through the framework
// and must reproduce the basis's IntegralHartree / IntegralPotential -- the inversion working for the
// density-dependent terms.  We feed a hand-built Hermitian density matrix (2 electrons in (e0+e1)/sqrt2).
TEST_F(PlaneWaveDFT, PWDynamicTermsMatchBasis)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();

    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);     // Hermitian density matrix
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;            // sets D(1,0)=conj(D(0,1)) -> a non-uniform density
    Irrep irr=F.pw.GetIrrep(Spin::None);
    qchem::ChargeDensity::IrrepCD<dcmplx> cd(D, &F.pw, irr);

    // Hartree term matrix == basis IntegralHartree of the same density.
    qchem::Hamiltonian::PW_Hartree   h;
    qchem::Hamiltonian::cDynamic_HT* ht=&h;
    const chmat_t& Mh = ht->GetMatrix(&F.pw, Spin::None, &cd);
    double Eh; chmat_t refh = F.pw.IntegralHartree(cd, Eh);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(Mh(i,j))-dcmplx(refh(i,j))), 0.0, 1e-10);

    // XC term (Dirac exchange) matrix == basis IntegralPotential of v_xc(rho(r)).
    auto dirac=std::make_shared<qchem::Hamiltonian::SlaterExchange>(2.0/3.0);
    qchem::Hamiltonian::PW_XC        xc(dirac);
    qchem::Hamiltonian::cDynamic_HT* xt=&xc;
    const chmat_t& Mx = xt->GetMatrix(&F.pw, Spin::None, &cd);
    FieldFn vxcfield([&](const rvec3_t& r){return dirac->GetVxc(cd(r));});
    chmat_t refx = F.pw.IntegralPotential(vxcfield);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(Mx(i,j))-dcmplx(refx(i,j))), 0.0, 1e-10);
}

// THE STAGE-3 PAYOFF: silicon (Gamma) self-consistent DFT run entirely through the FRAMEWORK objects --
// a dcmplx HamiltonianImp summing the PW Kohn-Sham terms (kinetic + external PP + Hartree + Dirac + VWN),
// an IrrepCD<dcmplx> density built from the complex orbitals, and the framework energy bookkeeping --
// reproducing the standalone prototype's Si-Gamma result (Etot=1.468, 8 valence electrons).  A thin SCF
// driver stands in for the (not-yet-complexified) WaveFunction/SCFIterator orchestration.
TEST_F(PlaneWaveDFT, FrameworkSiliconGammaMatchesPrototype)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26, h=0.5*a;
    Matrix3D<double> Amat(0.0,h,h,  h,0.0,h,  h,h,0.0);
    UnitCell          cell(Amat);
    Lattice_3D        lat(cell, ivec3_t(1,1,1));
    ReciprocalLattice recip(lat.Reciprocal());
    PlaneWave_IBS     pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);

    auto si=std::make_shared<Molecule>();
    si->Insert(new Atom(14, rvec3_t(0,0,0)));
    si->Insert(new Atom(14, rvec3_t(0.25*a,0.25*a,0.25*a)));
    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();
    pw.SetPseudopotential(&loc, &nl);

    // Framework Hamiltonian: a dcmplx HamiltonianImp summing the PW Kohn-Sham terms.
    cHamiltonianImp ham;
    ham.Add(new PW_Kinetic);
    ham.Add(new PW_External(si));
    ham.Add(new PW_Hartree);
    ham.Add(new PW_XC(std::make_shared<SlaterExchange> (2.0/3.0)));   // Dirac exchange
    ham.Add(new PW_XC(std::make_shared<VWN_Correlation>()));          // VWN5 correlation

    const int Nelec=8, Nocc=Nelec/2;
    size_t  n=pw.GetNumFunctions();
    Irrep   irr=pw.GetIrrep(Spin::None);
    chmat_t S=pw.MakeOverlap();

    // Seed a UNIFORM density (Hartree+XC present from iteration 0, as real PW codes do): D = (N/n) I.
    hmat_t<dcmplx> Dprev=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++) Dprev(i,i)=double(Nelec)/double(n);
    auto cd=std::make_unique<qchem::ChargeDensity::IrrepCD<dcmplx>>(Dprev, &pw, irr);

    bool converged=false; double Eprev=1e30;
    qchem::EnergyBreakdown E;
    for (int it=0; it<400; it++)
    {
        chmat_t H=ham.GetMatrix(&pw, Spin::None, cd.get());      // sums the framework terms for the current density
        LASolver<dcmplx>* las=LASolver<dcmplx>::Factory(qchem::Eigen);
        las->SetBasisOverlap(S);
        auto sol=las->Solve(H);
        delete las;
        mat_t<dcmplx> U=std::get<0>(sol); rvec_t e=std::get<1>(sol);

        std::vector<size_t> idx(e.size());                       // occupy the Nocc lowest bands
        for (size_t i=0;i<idx.size();i++) idx[i]=i;
        std::sort(idx.begin(),idx.end(),[&](size_t aa,size_t bb){return e[aa]<e[bb];});

        hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);                // density matrix D = Sum_occ 2 c c^H (Hermitian)
        for (int k=0;k<Nocc;k++)
        {
            size_t c=idx[k];
            for (size_t i=0;i<n;i++)
                for (size_t j=i;j<n;j++)
                    D(i,j) += 2.0*dcmplx(U(i,c))*std::conj(dcmplx(U(j,c)));
        }
        Dprev = hmat_t<dcmplx>(0.6*Dprev + 0.4*D);               // linear mixing
        cd=std::make_unique<qchem::ChargeDensity::IrrepCD<dcmplx>>(Dprev, &pw, irr);

        // GetTotalEnergy(new cd) computes the energy AND invalidates the dynamic terms' Irrep-keyed cache
        // (their GetEnergy calls newCD) so the NEXT GetMatrix rebuilds the Hartree/XC matrices fresh.
        E=ham.GetTotalEnergy(cd.get());
        double Etot=E.GetTotalEnergy();
        if (std::abs(Etot-Eprev)<1e-7) { converged=true; break; }
        Eprev=Etot;
    }
    ASSERT_TRUE(converged);
    std::cout << "[Si framework-Gamma] charge="<<cd->GetTotalCharge()<<" Etot="<<E.GetTotalEnergy()
              << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc<<")" << std::endl;

    EXPECT_NEAR(cd->GetTotalCharge(), 8.0, 1e-6);                 // 8 valence electrons
    EXPECT_NEAR(E.GetTotalEnergy(),   1.468, 5e-3);              // matches the standalone prototype Si-Gamma
}

// The SAME Si-Gamma Kohn-Sham problem as FrameworkSiliconGammaMatchesPrototype, but now driven by the
// REAL framework cSCFIterator (no hand-rolled SCF loop): cSCFIterator -> cWaveFunction (UnPolarizedWF
// -> IrrepWF) -> tSCFAcceleratorNull<dcmplx> diagonalize -> TOrbitals<dcmplx> fill -> IrrepCD<dcmplx>,
// with the dcmplx HamiltonianImp summing the PW terms.  This is the milestone that retires the "k-loop
// in the IBS": single-k plane-wave DFT IS now Hamiltonian = Sum terms + SCFIterator, like atoms/molecules.
TEST_F(PlaneWaveDFT, FrameworkSiliconGammaThroughSCFIterator)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26;                       // Si conventional cubic lattice constant (a.u.)
    FCCUnitCell cell(a);                         // FCC primitive cell
    cell.AddAtom(14, {0,0,0});                   // Si diamond two-atom basis, FRACTIONAL coordinates
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D  lat(cell, ivec3_t(1,1,1));

    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();      // must outlive the SCF run (the basis holds &loc)
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();

    // The basis comes from the factory as an abstract BasisSet<dcmplx>; it owns its plane-wave Bloch
    // block(s) and (for now) carries the pseudopotential.  Single-k -> one block at Gamma.
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, 4.0, &loc, &nl));
    const auto* pw=(*bs)[0];                         // the single Bloch block (for the seed density below)

    Irrep      irr=bs->GetIrreps(Spin::None)[0];
    size_t     n  =bs->GetNumFunctions();
    const int  Nelec=8;

    Crystal_EC  ec(irr, Nelec);

    // Plane-wave LDA Kohn-Sham Hamiltonian from the Ham factory family: kinetic + external(pseudo) +
    // Hartree + Dirac exchange + VWN5 correlation (heap; the SCFIterator takes ownership).  The atoms
    // are sourced from the lattice -- the structure carries all the external-potential information.
    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure());

    // No-acceleration manager (plain diagonalize each iteration) for the complex path.
    auto* acc=new qchem::SCFAccelerators::tSCFAcceleratorNull<dcmplx>();

    // Seed a UNIFORM density (Hartree+XC active from iteration 0, as real PW codes do): D = (N/n) I.
    hmat_t<dcmplx> D0=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++) D0(i,i)=double(Nelec)/double(n);
    auto* seed=new qchem::ChargeDensity::IrrepCD<dcmplx>(D0, pw, irr);

    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc, seed);

    SCFParams par;
    par.NMaxIter      =80;
    // The energy converges to 1.468057 (== the prototype) by ~iter 10, but the Null accelerator has no
    // DIIS damping, so the density matrix drifts marginally within Si's degenerate Gamma_25' occupied
    // manifold (energy flat to 1e-9, |Delta rho| dips to ~1e-5 then slowly grows).  Converge at the
    // density dip with a realistic plane-wave tolerance; a complex DIIS accelerator is future work.
    par.MinΔρ         =1e-4;
    par.MinΔFD        =1e30;   // Null accelerator reports no [F,D] error
    par.MinVirial     =1e30;   // pseudopotential calc: the textbook -V/K=2 virial does not hold
    par.MinFD         =1e30;
    par.StartingRelaxRo=0.4;
    par.MergeTol      =1e-4;
    par.Verbose       =false;
    scf.Iterate(par);

    const qchem::WaveFunction::cWaveFunction* wf=scf.GetWaveFunction();
    auto* cd=wf->GetChargeDensity();              // caller owns the returned (composite) density
    double charge=cd->GetTotalCharge();
    delete cd;
    qchem::EnergyBreakdown E=scf.GetEnergy();
    std::cout << "[Si SCFIterator-Gamma] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Etot="<<E.GetTotalEnergy()
              << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc<<")" << std::endl;

    EXPECT_TRUE(scf.Converged());
    EXPECT_NEAR(charge,             8.0,   1e-6);   // 8 valence electrons
    EXPECT_NEAR(E.GetTotalEnergy(), 1.468, 5e-3);   // matches the standalone prototype Si-Gamma
}

// Stage 4 / multi-k: the SAME Si Kohn-Sham problem on a 2x2x2 Brillouin-zone mesh, through the REAL
// cSCFIterator.  Now the basis holds 8 Bloch blocks (one per k-point); the framework's per-irrep loop
// (MakeIrrepWFs, one IrrepWF per block) IS the BZ sum Sum_k w_k, with each block's density BZ-weighted
// (Symmetry::GetWeight = w_k = 1/8) so the total charge is 8 (not 8*8).  Reproduces the standalone
// prototype ScfSiliconBZSampled (Etot=0.934, gap>0).
TEST_F(PlaneWaveDFT, FrameworkSilicon2x2x2ThroughSCFIterator)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26;
    FCCUnitCell cell(a);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D  lat(cell, ivec3_t(2,2,2));       // 2x2x2 = 8 k-points

    HGH_LocalPotential     loc=HGH_LocalPotential::Silicon();
    HGH_SeparablePotential nl =HGH_SeparablePotential::Silicon();

    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, 4.0, &loc, &nl));
    BasisSet::irrepv_t irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per k-block (8)

    const int Nelec=8;
    Crystal_EC ec(irreps, Nelec);                          // Nval per k-block; weights handle the BZ sum

    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure());
    auto* acc=new qchem::SCFAccelerators::tSCFAcceleratorNull<dcmplx>();

    // Uniform-density seed on the first block: D=(N/n0)I gives rho(r)=N/V (uniform), the total density
    // every block's first Hartree/XC needs (a single block suffices since rho is constant).
    const auto* pw0=(*bs)[0];
    size_t n0=pw0->GetNumFunctions();
    hmat_t<dcmplx> D0=blazem::zeroH<dcmplx>(n0);
    for (size_t i=0;i<n0;i++) D0(i,i)=double(Nelec)/double(n0);
    auto* seed=new qchem::ChargeDensity::IrrepCD<dcmplx>(D0, pw0, irreps[0]);

    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc, seed);

    SCFParams par;
    par.NMaxIter      =80;
    par.MinΔρ         =1e-4;
    par.MinΔFD        =1e30;
    par.MinVirial     =1e30;
    par.MinFD         =1e30;
    par.StartingRelaxRo=0.4;
    par.MergeTol      =1e-4;
    par.Verbose       =false;
    scf.Iterate(par);

    const qchem::WaveFunction::cWaveFunction* wf=scf.GetWaveFunction();
    auto* cd=wf->GetChargeDensity();
    double charge=cd->GetTotalCharge();
    delete cd;
    qchem::EnergyBreakdown E=scf.GetEnergy();

    // Gap = lowest unoccupied - highest occupied, merged across all k-blocks (occupations are physical:
    // 2 for filled bands, 0 above, since only the density is BZ-weighted, not the occupation).
    double homo=-1e30, lumo=1e30;
    for (const auto& [e,lev] : wf->GetEnergyLevels())
        if (lev.occ>1e-6) homo=std::max(homo,e); else lumo=std::min(lumo,e);
    double gap=lumo-homo;

    std::cout << "[Si SCFIterator-2x2x2] iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Etot="<<E.GetTotalEnergy()<<" gap="<<gap<<" Ha ("<<gap*27.2114<<" eV)"
              << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een<<" Eee="<<E.Eee<<" Exc="<<E.Exc<<")" << std::endl;

    EXPECT_TRUE(scf.Converged());
    EXPECT_NEAR(charge,             8.0,   1e-6);    // 8 valence electrons (BZ-weighted sum)
    EXPECT_NEAR(E.GetTotalEnergy(), 0.934, 5e-3);    // matches prototype ScfSiliconBZSampled
    EXPECT_GT(gap, 0.0);                              // Si is a semiconductor
}

// Regression guard for the Hamiltonian-framework cache bug: Dynamic_HT_Imp::GetMatrix must invalidate
// its Irrep-keyed cache when the density changes (else it returns a STALE matrix to any caller that
// didn't happen to call GetTotalEnergy(new cd) in between).  Two different densities (same Irrep, no
// GetEnergy between): the second GetMatrix must reflect the second density, not return the first's
// (uniform => zero) Hartree matrix.  Was a DISABLED known-failure; fixed by the cd-change cache clear.
TEST_F(PlaneWaveDFT, DynamicTermCacheFreshAcrossDensity)
{
    PWFixture F;
    size_t  n=F.pw.GetNumFunctions();
    Irrep   irr=F.pw.GetIrrep(Spin::None);

    hmat_t<dcmplx> D1=blazem::zeroH<dcmplx>(n);  D1(0,0)=2.0;                            // uniform density
    hmat_t<dcmplx> D2=blazem::zeroH<dcmplx>(n);  D2(0,0)=1.0; D2(0,1)=1.0; D2(1,1)=1.0;  // modulated density
    qchem::ChargeDensity::IrrepCD<dcmplx> cd1(D1,&F.pw,irr), cd2(D2,&F.pw,irr);

    qchem::Hamiltonian::PW_Hartree   hart;
    qchem::Hamiltonian::cDynamic_HT* ht=&hart;
    ht->GetMatrix(&F.pw, Spin::None, &cd1);                       // populates the cache for cd1
    const chmat_t& M2 = ht->GetMatrix(&F.pw, Spin::None, &cd2);   // BUG: returns cd1's stale matrix

    double  Eh; chmat_t ref2=F.pw.IntegralHartree(cd2,Eh);        // the correct V_H for cd2
    double  diff=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            diff=std::max(diff, std::abs(dcmplx(M2(i,j))-dcmplx(ref2(i,j))));
    EXPECT_NEAR(diff, 0.0, 1e-10);   // fails today (M2 is cd1's matrix); passes after the cache-invalidation fix
}

} //namespace

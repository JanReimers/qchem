// file: PlaneWaveDFTUT.C  Prototype self-consistent plane-wave DFT, validated outside the Hamiltonian
// framework (see doc/PlaneWavePlan.md sequencing + memory project_dft_upgrade_plan).
//
// WHY this lives in a unit test, not a library: the prototype wires three libraries together --
// BasisSet/Lattice_3D (PlaneWave_IBS primitives), Hamiltonian (the validated LDA ExFunctional), and
// LASolver (complex eigensolver).  A library that did this would invert the layering (BasisSet sits
// BELOW Hamiltonian).  The test is the one place allowed to reach across all three.  Once the
// Hamiltonian/SCFIterator stack is templated on T (double|dcmplx) these pieces migrate to their proper
// homes: the density -> a dcmplx ChargeDensity, Hartree/XC -> tDynamic_HT<dcmplx> terms.
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
#include <cstdio>
#include "gtest/gtest.h"

import qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.BasisSet.Lattice_3D.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
import qchem.Mesh;                           // qcMesh::MeshParams (PW_Hartree's fit-basis factory arg; ignored)
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
import qchem.Blaze;          // mat_t<dcmplx>
import qchem.Math;           // Pi
import qchem.Hamiltonian.Internal.ExFunctional;     // the validated LDA functional interface
import qchem.Hamiltonian.Internal.SlaterExchange;   // Dirac exchange (alpha=2/3), eps_x = 3/4 v_x
import qchem.Hamiltonian.Internal.VWN_Correlation;  // VWN5 correlation (validated vs libxc)
import qchem.Hamiltonian.Internal.PWTerms;          // PW_Pseudo (dcmplx Hamiltonian term)
import qchem.Hamiltonian;                           // cStatic_HT / cDynamic_HT aliases (public term interfaces)
import qchem.Hamiltonian.Internal.Hamiltonian;      // cHamiltonianImp (the dcmplx Hamiltonian = sum of terms)
import qchem.Hamiltonian.Internal.Hamiltonians;     // Ham_PW_DFT (the assembled plane-wave LDA KS Hamiltonian)
import qchem.Hamiltonian.Internal.Terms;            // Vnn (ion-ion term: pair sum / Ewald via isFinite)
import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
import qchem.Energy;                                // EnergyBreakdown
import qchem.ChargeDensity.Imp.IrrepCD;             // IrrepCD<dcmplx> (concrete complex density)
import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
import qchem.Math.GMap;                             // ΔG_Map (the G-space coefficient map)
import qchem.BasisSet.Fit_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
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
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;      // cSCFAcceleratorDIIS (complex DIIS)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)
using namespace qchem;

using BasisSet::Lattice_3D::PlaneWave_IBS;
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

namespace
{

// PW_Hartree/PW_XC now take their fit basis (from the basis's own factory) at construction, like
// FittedVee/FittedVxc.  These low-level term tests build one straight from the plane-wave basis at hand.
qchem::Hamiltonian::PW_Hartree* NewPWHartree(const PlaneWave_IBS& pw)
{
    return new qchem::Hamiltonian::PW_Hartree(
        qchem::Hamiltonian::PW_Hartree::fbs_t(pw.CreateCDFitBasisSet(nullptr, qcMesh::MeshParams{})));
}
qchem::Hamiltonian::PW_XC* NewPWXC(const PlaneWave_IBS& pw, const qchem::Hamiltonian::PW_XC::xc_t& xc)
{
    return new qchem::Hamiltonian::PW_XC(xc,
        qchem::Hamiltonian::PW_XC::fbs_t(pw.CreateVxcFitBasisSet(nullptr, qcMesh::MeshParams{})));
}

// A ScalarFunction<double> wrapping a lambda f(r) -- to hand real-space fields to the basis's
// high-level integral methods (Overlap/Repulsion/Integral).
struct FieldFn : public ScalarFunction<double>
{
    std::function<double(const rvec3_t&)> f;
    explicit FieldFn(std::function<double(const rvec3_t&)> g) : f(g) {}
    virtual double  operator()(const rvec3_t& r) const {return f(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const {return rvec3_t(0,0,0);}
};

// A v_xc(r) real field presented as a ProjectedScalar_R, so a scalar fitter can fit it (item K probe).
struct FieldFnR : public ScalarFunction<double>, public qchem::Fitting::ProjectedScalar_R
{
    std::function<double(const rvec3_t&)> f;
    explicit FieldFnR(std::function<double(const rvec3_t&)> g) : f(g) {}
    virtual double  operator()(const rvec3_t& r) const override {return f(r);}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const override {return rvec3_t(0,0,0);}
    virtual const ScalarFunction<double>* GetScalarFunction() const override {return this;}
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
// templated Hamiltonian/SCFIterator<dcmplx> will eventually replace (V_H/V_xc -> tDynamic_HT<dcmplx>).

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
                         const Structure& st, const HGH_LocalPotential& loc, const HGH_SeparablePotential& nl,
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
                Vext.push_back(p->MakeLocalPotential(&st,loc)+p->MakeSeparablePotential(&st,nl));
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

// Item K EXPLORATION on a SELF-CONSISTENT density: converge the weak-cosine LDA density, then sweep the
// Vxc FIT grid (relCutoff 1->4->16, i.e. 16^3->32^3->64^3) and print how the XC energy quadrature E_xc=∫ε ρ
// and the SPATIAL fit residual ‖v_xc - v_xc,fit‖ behave.  Illustrative (prints a table); the density is real.
TEST_F(PlaneWaveDFT, ItemK_Explore_ScfDensity)
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
    chmat_t Vext=pw.MakePotential([V0](const ivec3_t& dm)->dcmplx
    { return (dm.x*dm.x+dm.y*dm.y+dm.z*dm.z==1) ? dcmplx(V0) : dcmplx(0.0); });

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, 2, ivec3_t(12,12,12), vxcOf, epsOf, 0.5, 1e-9, 400);
    ASSERT_TRUE(R.converged);

    ΔG_Map rho;                                              // the self-consistent density as a ΔG_Map
    for (const auto& kv : R.rho) rho[kv.first]=kv.second;
    FieldFnR vtrue([&](const rvec3_t& r){return vxcOf(pw.EvalField(rho, r));});   // exact v_xc(r) of the SCF density

    std::vector<rvec3_t> ref;                                // fit-grid-INDEPENDENT reference points (8^3)
    for (int i=0;i<8;i++) for (int j=0;j<8;j++) for (int k=0;k<8;k++)
        ref.push_back(cell.ToCartesian(rvec3_t((i+0.5)/8.0,(j+0.5)/8.0,(k+0.5)/8.0)));

    std::cout << "\n  SCF: converged in " << R.iters << " iters, Etot_direct=" << R.Etot_direct
              << ", Exc(scf grid)=" << R.Exc << "\n"
              << "  relCutoff | nG_fit |   FFT pts |      E_xc(∫ε ρ) | ‖v_xc−v_xc,fit‖\n"
              <<   "  ----------+--------+-----------+-----------------+----------------\n";
    std::vector<double> Excs, resids;
    for (double rc : {1.0, 2.0, 4.0})
    {
        qcMesh::MeshParams mp; mp.relCutoff=rc;
        auto fb=qchem::Hamiltonian::PW_XC::fbs_t(pw.CreateVxcFitBasisSet(nullptr, mp));
        auto ge=dynamic_cast<const qchem::BasisSet::G_FieldEvaluator*>(fb.get());
        size_t nG=fb->GetNumFunctions(), Npts=ge->GridPoints().size();

        rvec_t rgrid=ge->RhoOnGrid(rho), exc(rgrid.size());
        for (size_t q=0;q<rgrid.size();q++) exc[q]=epsOf(rgrid[q])*rgrid[q];
        double Exc=ge->Integral(exc);

        auto fitter=qchem::Fitting::Factory(fb);
        fitter->DoFit(vtrue);
        const ScalarFunction<double>& vfit=*fitter;
        double s=0.0;
        for (const rvec3_t& r : ref){double d=vfit(r)-vtrue(r); s+=d*d;}
        double resid=std::sqrt(s/ref.size());

        char row[256];
        std::snprintf(row,sizeof row,"  %8.0f  | %6zu | %9zu | %15.10f | %14.6e\n",rc,nG,Npts,Exc,resid);
        std::cout << row;
        Excs.push_back(Exc); resids.push_back(resid);
    }
    std::cout << std::endl;

    // The lesson, asserted: the SPATIAL fit residual converges FAST (spectrally, for this smooth density) as
    // the grid densifies -- but the XC ENERGY is essentially BLIND to it (flat to <1e-6).  So acceptance MUST
    // be field/density convergence vs a fine reference, NEVER dE_total (the non-variational-fit pin).
    EXPECT_LT(resids[1], resids[0]);              // residual converges as the fit grid densifies ...
    EXPECT_LT(resids[2], resids[1]);              // ... monotonically
    EXPECT_LT(std::abs(Excs[2]-Excs[0]), 1e-6);   // while E_xc barely moves -- energy is the wrong yardstick
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

    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;
    chmat_t Vext = pw.MakeLocalPotential(&si,loc) + pw.MakeSeparablePotential(&si,nl);

    qchem::Hamiltonian::SlaterExchange  ex(2.0/3.0);
    qchem::Hamiltonian::VWN_Correlation vwn;
    auto vxcOf=[&](double r){return ex.GetVxc (r)+vwn.GetVxc (r);};
    auto epsOf=[&](double r){return ex.GetEpsXc(r)+vwn.GetEpsXc(r);};

    const int Nval=8;                                      // 2 Si x Zion 4
    int       m=MaxGComponent(pw);
    ivec3_t   Ng(4*m+1, 4*m+1, 4*m+1);                     // resolves the difference set (no aliasing)

    SCFResult R=RunSCF(pw, recip.GetCell(), Omega, Vext, Nval, Ng, vxcOf, epsOf, 0.4, 1e-8, 400);

    // Ion-ion (Ewald) energy of the two Si4+ cores: turns the electronic energy (which drops the G=0
    // potential and has no ion-ion term) into a physical, NEGATIVE total.  The Ewald cell carries the
    // same FCC geometry and the diamond basis given in fractional coordinates; charges = Zion = 4.
    UnitCell ecell(A);
    ecell.AddAtom(14, rvec3_t(0.0,0.0,0.0));
    ecell.AddAtom(14, rvec3_t(0.25,0.25,0.25));
    double Eii  = EwaldEnergy(ecell, rvec_t{4.0,4.0});

    // G=0 alignment of the local pseudopotential: the dropped G=0 potential carries a finite
    // electron-ion shift E_alpha = (N/Omega) Sum_a alpha_a, alpha_a = integral[V_loc^a + Zion/r].
    // (Two Si atoms, same species.)  This is the last G=0 piece needed for the absolute total energy.
    double Ealpha = (double(Nval)/Omega) * 2.0 * loc.FormFactorG0(14);
    double Etot   = R.Etot_direct + Eii + Ealpha;

    std::cout << "[Si] nG="<<pw.GetNumFunctions()<<" grid="<<(4*m+1)<<"^3 iters="<<R.iters
              << " converged="<<R.converged << "\n  Ekin="<<R.Ekin<<" Eext="<<R.Eext
              << " E_H="<<R.EH<<" E_xc="<<R.Exc << "\n  Etot(elec)="<<R.Etot_direct
              << " E_ion-ion(Ewald)="<<Eii << " E_alpha(G=0)="<<Ealpha
              << "\n  Etot(physical)="<<Etot << std::endl;

    ASSERT_TRUE(R.converged);
    EXPECT_NEAR(Omega*std::real(RhoAt(R.rho,ivec3_t(0,0,0))), double(Nval), 1e-6);  // 8 valence e-
    EXPECT_NEAR(R.Etot_band, R.Etot_direct, 1e-5);                                  // stationary fixed point
    EXPECT_LT(Etot, 0.0);   // ion-ion Madelung + G=0 alignment make the total energy negative
}

// The generalized Vnn Hamiltonian term: for a periodic Structure (isFinite()==false) it routes the
// ion-ion energy through the Ewald lattice sum (charges = the cell atoms' itsZ = the ion/valence charge),
// rather than the conditionally-convergent direct pair sum used for finite molecules.
TEST_F(PlaneWaveDFT, VnnPeriodicUsesEwald)
{
    const double a=10.26, h=0.5*a;
    Matrix3D<double> A(0.0,h,h,  h,0.0,h,  h,h,0.0);   // FCC primitive
    auto cell=std::make_shared<UnitCell>(A);
    cell->AddAtom(4,rvec3_t(0.0,0.0,0.0));             // itsZ = Zion = 4 (Si valence)
    cell->AddAtom(4,rvec3_t(0.25,0.25,0.25));
    std::shared_ptr<const Structure> st=cell;
    EXPECT_FALSE(st->isFinite());

    qchem::Hamiltonian::Vnn vnn(st);
    EnergyBreakdown eb;
    vnn.GetEnergy(eb, nullptr);                        // periodic branch ignores the density
    double ref=EwaldEnergy(*cell, rvec_t{4.0,4.0});
    EXPECT_NEAR(eb.Enn, ref, 1e-9);                    // routes through Ewald
    EXPECT_NEAR(eb.Enn, -8.40046, 1e-4);              // == the Si ion-ion Madelung energy
}

// (The GTH database-reader unit test lives in src/BasisSet/Lattice_3D/tests/GTH_UT.C -- it needs only
// the basis layer, not the SCF stack.)

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
    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;

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

// --- Stage 2: the basis-level high-level DFT capability (Band_DFT_IBS), validated against the
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

// Repulsion on a single-cosine density rho = rho0 + 2A cos(G0.r) reproduces the Poisson result
// (E_H = Omega 4 pi A^2/|G0|^2, V_H~(G0) = 4 pi A/|G0|^2) -- same check as HartreeSingleCosineMatchesPoisson,
// now driven entirely through the basis's high-level method (the basis samples + transforms internally).
TEST_F(PlaneWaveDFT, BasisRepulsionMatchesPoisson)
{
    PWFixture F;
    const double rho0=2.0/F.Omega, A=0.01;
    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    double  g02=G0*G0;
    FieldFn rho([&](const rvec3_t& r){return rho0 + 2*A*std::cos(G0*r);});

    double  Eh=0.0;
    chmat_t VH=F.pw.Repulsion(rho, Eh);
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

// Overlap of a constant field f is f*Identity; of a reciprocal cosine has f-tilde(+/-e)=1/2.
TEST_F(PlaneWaveDFT, BasisIntegralPotentialConstAndCosine)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();

    FieldFn cst([](const rvec3_t&){return 0.7;});
    chmat_t Vc=F.pw.Overlap(cst);
    for (size_t i=0;i<n;i++)
    {
        EXPECT_NEAR(std::real(dcmplx(Vc(i,i))), 0.7, 1e-9);
        for (size_t j=i+1;j<n;j++) EXPECT_NEAR(std::abs(dcmplx(Vc(i,j))), 0.0, 1e-9);
    }

    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    FieldFn cosfn([&](const rvec3_t& r){return std::cos(G0*r);});
    chmat_t Vk=F.pw.Overlap(cosfn);
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

// The PW_Pseudo Hamiltonian term (cStatic_HT<dcmplx>) routes through the framework's GetMatrix path
// and the abstract Band_DFT_IBS dynamic_cast -- its matrix must equal the basis's external assembly.
// This is the dependency inversion working end-to-end through a real Hamiltonian term.
TEST_F(PlaneWaveDFT, PWPseudoTermMatchesBasis)
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
    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;

    // The TERM owns the model now and asks the basis to assemble it (the pseudo-wall).
    qchem::Hamiltonian::PW_Pseudo   ext(si, &loc, &nl);
    qchem::Hamiltonian::cStatic_HT*   term=&ext;        // the public term interface (as the Hamiltonian holds it)
    const chmat_t& M  = term->GetMatrix(&pw, Spin::None);
    chmat_t        ref= pw.MakeLocalPotential(si.get(),loc) + pw.MakeSeparablePotential(si.get(),nl);

    size_t n=pw.GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(M(i,j))-dcmplx(ref(i,j))), 0.0, 1e-12);
}

// The PW_Hartree and PW_XC dynamic terms route a complex density (IrrepCD<dcmplx>) through the framework
// and must reproduce the basis's Repulsion / Overlap -- the inversion working for the
// density-dependent terms.  We feed a hand-built Hermitian density matrix (2 electrons in (e0+e1)/sqrt2).
TEST_F(PlaneWaveDFT, PWDynamicTermsMatchBasis)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();

    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);     // Hermitian density matrix
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;            // sets D(1,0)=conj(D(0,1)) -> a non-uniform density
    Irrep irr=F.pw.GetIrrep(Spin::None);
    qchem::ChargeDensity::IrrepCD<dcmplx> cd(D, &F.pw, irr);

    // Hartree term matrix == basis Repulsion of the same density.
    std::unique_ptr<qchem::Hamiltonian::PW_Hartree> h(NewPWHartree(F.pw));
    qchem::Hamiltonian::cDynamic_HT* ht=h.get();
    const chmat_t& Mh = ht->GetMatrix(&F.pw, Spin::None, &cd);
    double Eh; chmat_t refh = F.pw.Repulsion(cd, Eh);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(Mh(i,j))-dcmplx(refh(i,j))), 0.0, 1e-10);

    // XC term (Dirac exchange) matrix == the basis's FFT route: rho(r) via inverse FFT of rho-tilde,
    // v_xc applied pointwise on the grid, forward FFT to the matrix (what the term itself does).
    auto dirac=std::make_shared<qchem::Hamiltonian::SlaterExchange>(2.0/3.0);
    std::unique_ptr<qchem::Hamiltonian::PW_XC> xc(NewPWXC(F.pw, dirac));
    qchem::Hamiltonian::cDynamic_HT* xt=xc.get();
    const chmat_t& Mx = xt->GetMatrix(&F.pw, Spin::None, &cd);
    rvec_t rho=F.pw.RhoOnGrid(F.pw.MakeFourierDensity(D));
    rvec_t vxc(rho.size());
    for (size_t q=0;q<rho.size();q++) vxc[q]=dirac->GetVxc(rho[q]);
    chmat_t refx = F.pw.Overlap(vxc);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            EXPECT_NEAR(std::abs(dcmplx(Mx(i,j))-dcmplx(refx(i,j))), 0.0, 1e-10);
}

// Item K: the XC quadrature grid now comes from the FIT basis, so relCutoff (the CP2K REL_CUTOFF density-
// cutoff idea) is the genuine fit-accuracy lever.  (1) a larger relCutoff yields a DENSER fit {G}; (2) it
// feeds through to the Vxc matrix (relCutoff=1 -> the orbital grid, bit-identical to the FFT route; >1 -> a
// finer, less-aliased quadrature); (3) the matrix GRID-CONVERGES (successive refinements shrink).  The fit is
// non-variational, so acceptance is convergence of the matrix, NEVER an energy anchor.
TEST_F(PlaneWaveDFT, ItemK_RelCutoffDensifiesAndConvergesVxc)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;                       // non-uniform => nonlinear v_xc aliases on a coarse grid
    Irrep irr=F.pw.GetIrrep(Spin::None);
    qchem::ChargeDensity::IrrepCD<dcmplx> cd(D, &F.pw, irr);
    auto dirac=std::make_shared<qchem::Hamiltonian::SlaterExchange>(2.0/3.0);

    auto vxcAt=[&](double relCutoff, size_t& nGfit)
    {
        qcMesh::MeshParams mp; mp.relCutoff=relCutoff;
        auto fb=qchem::Hamiltonian::PW_XC::fbs_t(F.pw.CreateVxcFitBasisSet(nullptr, mp));
        nGfit=fb->GetNumFunctions();
        std::unique_ptr<qchem::Hamiltonian::PW_XC> xc(new qchem::Hamiltonian::PW_XC(dirac, fb));
        return chmat_t(static_cast<qchem::Hamiltonian::cDynamic_HT*>(xc.get())->GetMatrix(&F.pw, Spin::None, &cd));
    };
    auto froDiff=[&](const chmat_t& A, const chmat_t& B)
    {
        double s=0.0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) s+=std::norm(dcmplx(A(i,j))-dcmplx(B(i,j)));
        return std::sqrt(s);
    };

    size_t n1,n4,n16;
    chmat_t M1=vxcAt(1.0,n1), M4=vxcAt(4.0,n4), M16=vxcAt(16.0,n16);

    // (1) relCutoff densifies the fit {G}.
    EXPECT_GT(n4, n1);
    EXPECT_GT(n16, n4);

    // (2) relCutoff=1 reproduces the orbital-grid FFT route exactly (the fix is inert at Gamma).
    rvec_t rho=F.pw.RhoOnGrid(F.pw.MakeFourierDensity(D));
    rvec_t vxc(rho.size());
    for (size_t q=0;q<rho.size();q++) vxc[q]=dirac->GetVxc(rho[q]);
    chmat_t refx=F.pw.Overlap(vxc);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++)
        EXPECT_NEAR(std::abs(dcmplx(M1(i,j))-dcmplx(refx(i,j))), 0.0, 1e-10);

    // (3) the denser grid CHANGES the matrix (relCutoff is live) and it GRID-CONVERGES (Cauchy: the later,
    //     finer refinement moves the matrix less than the earlier one).
    double d1=froDiff(M4,M1), d2=froDiff(M16,M4);
    EXPECT_GT(d1, 1e-9);      // the fit grid actually feeds through to the quadrature
    EXPECT_LT(d2, d1);        // convergence
}

// Item B: the ortho (plane-wave) density fitter's fitted field is a REAL, evaluatable ScalarFunction --
// rho_fit(r) = Re Σ_dm rho~(dm) e^{i(B·dm)·r} (what the GUI plots).  Fit a non-uniform density through the
// production Factory path, then check the fitter's op(r)/Gradient via the G_FieldEvaluator DIP seam against
// (a) an independent inverse transform, (b) the FFT route RhoOnGrid at r=0, (c) finite-difference gradient.
TEST_F(PlaneWaveDFT, OrthoFitterRealSpaceField)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;                   // Hermitian, non-uniform -> a non-trivial rho~(dm)
    ΔG_Map rhoTilde=F.pw.MakeFourierDensity(D);

    // The Factory-built ortho density fitter (the production path); the core IS-A ScalarFunction<double> now.
    auto fb=std::shared_ptr<const qchem::BasisSet::cFIT_CD_ABS>(F.pw.CreateCDFitBasisSet(nullptr, qcMesh::MeshParams{}));
    auto fitter=qchem::Fitting::Factory(fb);
    fitter->DoFit(qchem::Fitting::ProjectedDensity_G(rhoTilde));
    const ScalarFunction<double>& rhoFit=*fitter;         // item B: FunctionFitter_Density<dcmplx> : ScalarFunction

    const UnitCell& B=F.B();
    auto direct=[&](const rvec3_t& r)                     // independent inverse transform (B.ToCartesian, not GetGCartesian)
    {
        dcmplx s(0.0);
        for (const auto& kv:rhoTilde){rvec3_t G=B.ToCartesian(rvec3_t(kv.first)); double ph=G*r; s+=kv.second*dcmplx(std::cos(ph),std::sin(ph));}
        return s.real();
    };
    for (const rvec3_t& r : {rvec3_t(0.0,0.0,0.0), rvec3_t(1.0,2.0,0.5), rvec3_t(3.1,1.7,2.4)})
        EXPECT_NEAR(rhoFit(r), direct(r), 1e-11);

    // Cross-check against the independent FFT route: RhoOnGrid[0] = rho(r=0) = Σ rho~(dm).
    rvec_t grid=F.pw.RhoOnGrid(rhoTilde);
    EXPECT_NEAR(rhoFit(rvec3_t(0.0,0.0,0.0)), grid[0], 1e-11);

    // Gradient: finite-difference the field at a generic point.
    rvec3_t r0(1.3,0.7,2.1), g=rhoFit.Gradient(r0);
    const double h=1e-5;
    EXPECT_NEAR(g.x, (rhoFit(rvec3_t(r0.x+h,r0.y,r0.z))-rhoFit(rvec3_t(r0.x-h,r0.y,r0.z)))/(2*h), 1e-4);
    EXPECT_NEAR(g.y, (rhoFit(rvec3_t(r0.x,r0.y+h,r0.z))-rhoFit(rvec3_t(r0.x,r0.y-h,r0.z)))/(2*h), 1e-4);
    EXPECT_NEAR(g.z, (rhoFit(rvec3_t(r0.x,r0.y,r0.z+h))-rhoFit(rvec3_t(r0.x,r0.y,r0.z-h)))/(2*h), 1e-4);
}

// THE STAGE-3 PAYOFF: silicon (Gamma) self-consistent DFT run entirely through the FRAMEWORK objects --
// a cHamiltonianImp summing the PW Kohn-Sham terms (kinetic + external PP + Hartree + Dirac + VWN),
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
    GTH_PP                 siPP=GetGTH("Si","LDA",4);          // CP2K GTH-LDA q4, from the database
    const HGH_LocalPotential&     loc=siPP.local;
    const HGH_SeparablePotential& nl =siPP.nonlocal;

    // Framework Hamiltonian: a cHamiltonianImp summing the PW Kohn-Sham terms.  The external term
    // owns the pseudopotential model (the pseudo-wall) and assembles it through the basis.
    cHamiltonianImp ham;
    ham.Add(new PW_Kinetic);
    ham.Add(new PW_Pseudo(si, &loc, &nl));
    ham.Add(NewPWHartree(pw));
    ham.Add(NewPWXC(pw, std::make_shared<SlaterExchange> (2.0/3.0)));   // Dirac exchange
    ham.Add(NewPWXC(pw, std::make_shared<VWN_Correlation>()));          // VWN5 correlation

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
    // This manual Hamiltonian has no ion-ion term; the band-structure (electronic) energy is the
    // prototype anchor (the dropped-G=0 alignment Ealign is now separated out of it).
    EXPECT_NEAR(E.GetElectronicEnergy(), 1.468, 5e-3);           // matches the standalone prototype Si-Gamma
}

// The SAME Si-Gamma Kohn-Sham problem as FrameworkSiliconGammaMatchesPrototype, but now driven by the
// REAL framework cSCFIterator (no hand-rolled SCF loop): cSCFIterator -> cWaveFunction (UnPolarizedWF
// -> IrrepWF) -> tSCFAcceleratorNull<dcmplx> diagonalize -> TOrbitals<dcmplx> fill -> IrrepCD<dcmplx>,
// with the cHamiltonianImp summing the PW terms.  This is the milestone that retires the "k-loop
// in the IBS": single-k plane-wave DFT IS now Hamiltonian = Sum terms + SCFIterator, like atoms/molecules.
TEST_F(PlaneWaveDFT, FrameworkSiliconGammaThroughSCFIterator)
{
    using namespace qchem::Hamiltonian;
    const double a=10.26;                       // Si conventional cubic lattice constant (a.u.)
    FCCUnitCell cell(a);                         // FCC primitive cell
    cell.AddAtom(14, {0,0,0});                   // Si diamond: itsZ = TRUE species Z = 14; Zion comes from the PP (ZionFn)
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D  lat(cell, ivec3_t(1,1,1));

    // Run the Si-Gamma SCF at plane-wave cutoff \a Ecut under a given SEED strategy, EVERYTHING ELSE
    // identical (config, Hamiltonian, complex-DIIS accelerator, convergence params) -- so the only difference
    // is the seed and the iteration counts are a fair head-to-head.  Returns iters/convergence/charge/energy.
    namespace L3=BasisSet::Lattice_3D;
    struct Run { size_t iters; bool conv; double charge; qchem::EnergyBreakdown E; };
    auto run=[&](double Ecut, qchem::ChargeDensity::SeedStrategy seed)
    {
        using namespace qchem::Hamiltonian;
        // The basis is an abstract tBasisSet<dcmplx> owning its plane-wave Bloch block(s); the PP model lives
        // on the Hamiltonian term (the pseudo-wall), NOT the basis.
        std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, Ecut));
        Irrep      irr=bs->GetIrreps(Spin::None)[0];
        Crystal_EC ec(irr, 8);   // Si insulator: 8 valence electrons (2 atoms x Zion 4)
        // Plane-wave LDA Kohn-Sham Hamiltonian: kinetic + external(pseudo) + Hartree + Dirac X + VWN5
        // (heap; the SCFIterator takes ownership).  Atoms from the lattice; the PP model lives on the term.
        cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), "Si", "LDA", 4);   // one-call: looks up + owns the GTH PP
        // Complex DIIS (Pulay): extrapolate the Fock matrix from the [F,D] history; damps the marginal drift
        // within Si's degenerate Gamma_25' manifold so it converges to a tight |Delta rho|.
        using qchem::SCFAccelerators::DIISParams;
        auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(DIISParams{8, 0.5, 1e-10, 1e-9});
        qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc, seed, lat.GetStructure().get());

        SCFParams par;
        par.NMaxIter      =80;
        par.MinΔρ         =1e-7;   // tight convergence (complex DIIS damps the degenerate-manifold drift)
        par.MinΔFD        =1e30;
        par.MinVirial     =1e30;   // pseudopotential calc: the textbook -V/K=2 virial does not hold
        par.MinFD         =1e30;
        par.StartingRelaxRo=0.4;
        par.MergeTol      =1e-4;
        par.Verbose       =false;
        scf.Iterate(par);

        auto* cd=scf.GetWaveFunction()->GetChargeDensity();   // caller owns the returned (composite) density
        double charge=cd->GetTotalCharge();
        delete cd;
        return Run{scf.GetIterationCount(), scf.Converged(), charge, scf.GetEnergy()};
    };

    // The plane-wave SAD path (FourierSeedCD: a G-space form-factor sum of the atomic VALENCE densities from
    // atomic_valence_densities.json).  That file now holds the SMOOTH pseudo-valence Si density produced by
    // the Atom-PP (Ham_PP + KB nonlocal: scfrun --model PP --valence 4 --out ...), not the old all-electron
    // valence whose core peak injected spurious high-G content.  See doc/SCFSeedingPlan.md section 9.7.

    // Ecut=4 / Gamma (fast, near-jellium).  Both seeds; cross-check the SAD energy vs the standalone prototype.
    Run sad=run(4.0, qchem::ChargeDensity::SeedStrategy::SAD);
    Run uni=run(4.0, qchem::ChargeDensity::SeedStrategy::Uniform);
    std::cout << "[Si SCFIterator-Gamma] SAD iters="<<sad.iters<<"  Uniform iters="<<uni.iters
              << "  SAD Etot="<<sad.E.GetTotalEnergy()
              << "  (Ekin="<<sad.E.Kinetic<<" Een="<<sad.E.Een<<" Eee="<<sad.E.Eee<<" Exc="<<sad.E.Exc<<")" << std::endl;

    EXPECT_TRUE(sad.conv);
    EXPECT_NEAR(sad.charge,                  8.0,    1e-6);   // 8 valence electrons
    EXPECT_NEAR(sad.E.GetElectronicEnergy(), 1.468,  5e-3);   // band energy matches the standalone prototype
    // Physical total: electronic + ion-ion Ewald (Enn) + dropped-G=0 alignment (Ealign) -- now NEGATIVE.
    // (Underconverged at Ecut=4 / Gamma-only; the converged Si total is ~-7.9 Ha/cell.)
    EXPECT_NEAR(sad.E.GetTotalEnergy(),     -7.2273, 5e-3) << "Enn="<<sad.E.Enn<<" Ealign="<<sad.E.Ealign;
    EXPECT_TRUE(uni.conv);
    EXPECT_NEAR(sad.E.GetTotalEnergy(), uni.E.GetTotalEnergy(), 1e-3);   // seed cannot change the answer

    // THE SMOOTH-DENSITY WIN (regression guard).  The smooth pseudo-valence seed converges in ~12 iters here;
    // the OLD all-electron-valence seed needed ~15 (its core peak injected spurious high-G content).  Guard
    // the robust 3-iter improvement -- NOT a brittle tie with Uniform: at this coarse near-jellium Ecut=4 the
    // converged density is nearly flat so Uniform (the pure G=0 seed) is already near-optimal (11).  The
    // STRICT win over Uniform appears once the density has real structure -- e.g. Ecut=9/Gamma gives SAD 10
    // vs Uniform 11 (verified manually) -- but a 1-iter margin behind an ~18s high-Ecut SCF is too slow and
    // brittle to gate CI.  See doc/SCFSeedingPlan.md section 9.7.
    EXPECT_LE(sad.iters, 13u) << "smooth pseudo-valence SAD should converge in <=13 iters (was ~15 all-electron)";
}

// Shared Gamma-point framework-SCF driver for the multi-species crystal tests: build the PW basis from
// the lattice, seed a UNIFORM density, run the full SCFIterator (complex DIIS), and report.  Takes
// ownership of \a ham (the SCFIterator deletes it).  Mirrors the Si framework test body.
namespace {
struct FwResult { bool converged; double charge; qchem::EnergyBreakdown E; size_t iters; };
// Takes a Ham FACTORY (not a pre-built Ham) so the Hamiltonian is built WITH the basis -- the fit-basis
// seam: Ham_PW_DFT now needs the composite basis to create the Hartree density-fit basis (like the
// molecular Ham DFT ctors take bs for FittedVee).
FwResult RunFrameworkGamma(const Lattice_3D& lat, double Ecut, int Nelec,
                           std::function<qchem::Hamiltonian::cHamiltonian*(const BasisSet::Complex_BS*)> mkHam,
                           const char* label,
                           qchem::ChargeDensity::SeedStrategy seed=qchem::ChargeDensity::SeedStrategy::Uniform)
{
    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, Ecut));
    qchem::Hamiltonian::cHamiltonian* ham=mkHam(bs.get());   // build the Ham WITH the basis (fit-basis seam)
    Irrep  irr=bs->GetIrreps(Spin::None)[0];
    size_t n  =bs->GetNumFunctions();
    Crystal_EC ec(irr, Nelec);
    using qchem::SCFAccelerators::DIISParams;
    // EMax must exceed the ionic [F,D] error (~1.4-3) or DIIS bails ("En>EMax") and linear mixing
    // oscillates on the strong Madelung field.  Engage DIIS from the start (EMax large).
    auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(DIISParams{10, 8.0, 1e-10, 1e-9});
    // \a seed defaults to Uniform rho(r)=N/V; IonicSAD pre-bakes the formal-charge transfer (Na+ + F-).
    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc, seed, lat.GetStructure().get());
    SCFParams par;
    par.NMaxIter=120; par.MinΔρ=1e-6; par.MinΔFD=1e30; par.MinVirial=1e30; par.MinFD=1e30;
    par.StartingRelaxRo=0.3; par.MergeTol=1e-4; par.Verbose=false;
    scf.Iterate(par);
    const qchem::WaveFunction::cWaveFunction* wf=scf.GetWaveFunction();
    auto* cd=wf->GetChargeDensity();
    double charge=cd->GetTotalCharge();
    delete cd;
    qchem::EnergyBreakdown E=scf.GetEnergy();
    std::cout << "["<<label<<"] nPW="<<n<<" iters="<<scf.GetIterationCount()<<" charge="<<charge
              << " Etot="<<E.GetTotalEnergy() << "  (Ekin="<<E.Kinetic<<" Een="<<E.Een
              << " Eee="<<E.Eee<<" Exc="<<E.Exc<<" Enn="<<E.Enn<<" Ealign="<<E.Ealign<<")" << std::endl;
    return {scf.Converged(), charge, E, scf.GetIterationCount()};
}
} // namespace

// Multi-species ionic crystal NaF (rocksalt = FCC + 2-atom basis), through the full SCFIterator with the
// multi-species Ham_PW_DFT facade.  Na (Zion=1) + F (Zion=7) = 8 valence electrons; the per-Z router model
// dispatches Na's vs F's pseudopotential per atom.  F's tight 2p sets the (high) cutoff.
TEST_F(PlaneWaveDFT, FrameworkNaFThroughSCFIterator)
{
    using namespace qchem::Hamiltonian;
    const double a=8.73;                          // NaF lattice constant ~4.62 A (a.u.)
    FCCUnitCell cell(a);                          // rocksalt: FCC lattice, 2-atom basis
    cell.AddAtom(11, {0,0,0});                    // Na  (true species Z=11; Zion=1 via the PP)
    cell.AddAtom(9,  {0.5,0.5,0.5});              // F   (true species Z=9;  Zion=7)
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    using SS=qchem::ChargeDensity::SeedStrategy;
    auto mkHam=[&](const BasisSet::Complex_BS* bs)   // multi-species one-call; bs -> the Hartree fit basis
        { return new Ham_PW_DFT(lat.GetStructure(), bs, {{"Na",1},{"F",7}}, "LDA"); };

    // IonicSAD seed: the electronegativity heuristic assigns Na+ / F-, and the per-species valence densities
    // are scaled to the formal-charge electron counts (0 e- on Na, 8 e- on F) -- the seed pre-bakes the
    // Na->F charge transfer and CONSERVES the cell electron count.  This asserts the seed is CORRECT (charge
    // sum rule + same converged energy as Uniform: the seed cannot change the answer).
    FwResult I=RunFrameworkGamma(lat, /*Ecut*/6.0, /*Nelec*/8, mkHam, "NaF IonicSAD", SS::IonicSAD);
    FwResult U=RunFrameworkGamma(lat, /*Ecut*/6.0, /*Nelec*/8, mkHam, "NaF Uniform",  SS::Uniform);
    std::cout << "[NaF seed compare] IonicSAD iters="<<I.iters<<"  Uniform iters="<<U.iters << std::endl;

    EXPECT_TRUE(I.converged);
    EXPECT_NEAR(I.charge, 8.0, 1e-6);             // 1 (Na) + 7 (F) valence electrons, conserved by the ionic seed
    // Regression anchor (Ecut=6, Gamma-only -> underconverged but deterministic), like the Si total.
    // Negative: the ionic Madelung (Enn~-14) + G=0 alignment dominate.
    EXPECT_NEAR(I.E.GetTotalEnergy(), -20.3293, 5e-3) << "Enn="<<I.E.Enn<<" Ealign="<<I.E.Ealign;
    EXPECT_NEAR(I.E.GetTotalEnergy(), U.E.GetTotalEnergy(), 1e-3);   // seed-independence of the converged answer
    // NOTE: IonicSAD does NOT (yet) beat Uniform on iterations for NaF (~69 vs ~58 at Ecut=6).  The crude
    // ionic density -- the neutral F valence scaled x8/7 -- is too COMPACT (the real F- electron is diffuse),
    // injecting high-G content much like the Si all-electron core peak did before the smooth pseudo-valence
    // density landed.  The iteration win awaits a proper anion (F-) pseudo-valence density; the heuristic +
    // charge-conserving seed wiring are correct and in place.  See doc/SCFSeedingPlan.md section 3.3.
}

// Heavy multi-species ionic crystal CsI (CsCl = simple cubic + 2-atom basis).  THE d-PROJECTOR TEST: both
// Cs (q1) and I (q7) carry l=2 (d) Kleinman-Bylander channels, the first DFT exercise of the (2l+1)P_2
// angular path (Si/NaF are l<=1).  Cs q1 (Zion=1) deliberately avoids the semicore q9, whose l=3 (f)
// projector the analytic HGH Qli table doesn't tabulate.  And the PP promise: Cs/I are SOFTER than F
// (bigger r_loc), so this heavy salt needs a LOWER cutoff than NaF.
TEST_F(PlaneWaveDFT, FrameworkCsIThroughSCFIterator)
{
    using namespace qchem::Hamiltonian;
    const double a=8.63;                          // CsI lattice constant ~4.567 A (a.u.)
    UnitCell cell(a);                             // CsCl: simple cubic, 2-atom basis
    cell.AddAtom(55, {0,0,0});                    // Cs  (true species Z=55; Zion=1, q1 -- s,p,d channels)
    cell.AddAtom(53, {0.5,0.5,0.5});              // I   (true species Z=53; Zion=7 -- s,p,d channels)
    Lattice_3D lat(cell, ivec3_t(1,1,1));

    auto mkHam=[&](const BasisSet::Complex_BS* bs)   // bs -> the Hartree fit basis (created once in BuildTerms)
        { return new Ham_PW_DFT(lat.GetStructure(), bs, {{"Cs",1},{"I",7}}, "LDA"); };

    FwResult R=RunFrameworkGamma(lat, /*Ecut*/4.0, /*Nelec*/8, mkHam, "CsI SCFIterator-Gamma");

    EXPECT_TRUE(R.converged);
    EXPECT_NEAR(R.charge, 8.0, 1e-6);             // 1 (Cs) + 7 (I) valence electrons (d-projectors active)
    // Regression anchor (Ecut=4, Gamma-only).  The point is the d-channel assembly runs end-to-end;
    // both species' l=2 Kleinman-Bylander projectors contribute via the (2l+1)P_2(cos gamma) path.
    EXPECT_NEAR(R.E.GetTotalEnergy(), -11.3868, 5e-3) << "Enn="<<R.E.Enn<<" Ealign="<<R.E.Ealign;
}

// G-space Hartree path: V_H assembled directly from the density's Fourier coefficients rho-tilde(dm)
// (PlaneWave_IBS::MakeFourierDensity -> Repulsion(ΔG_Map)) must equal the real-space route
// (sample rho(r) on the grid -> ForwardDFT -> V_H).  This is the O(n^2) G-space replacement for the
// O(Npts*n^2) pointwise sampling -- the foundation of the FFT speed-up.  Fast: no SCF, one block.
TEST_F(PlaneWaveDFT, HartreeFromFourierMatchesPointwise)
{
    const double a=10.26;
    FCCUnitCell    cell(a);
    Lattice_3D     lat(cell, ivec3_t(1,1,1));
    PlaneWave_IBS  pw(lat.Reciprocal(), ivec3_t(1,1,1), ivec3_t(0,0,0), 4.0);
    size_t n=pw.GetNumFunctions();
    ASSERT_GT(n, 3u);
    Irrep irr=pw.GetIrrep(Spin::None);

    // A non-trivial Hermitian density matrix: a uniform base plus some off-diagonal (cross-G) structure.
    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);
    for (size_t i=0;i<n;i++) D(i,i)=8.0/double(n);
    D(0,1)=dcmplx(0.10,0.05);
    D(0,2)=dcmplx(-0.07,0.02);
    D(1,2)=dcmplx(0.04,-0.06);
    qchem::ChargeDensity::IrrepCD<dcmplx> cd(D, &pw, irr);   // IS-A ScalarFunction rho(r)=phi^H D phi

    double EhA, EhB;
    chmat_t VA=pw.Repulsion(cd, EhA);                             // real-space: sample + ForwardDFT
    chmat_t VB=pw.Repulsion(pw.MakeFourierDensity(D), EhB);  // G-space: direct from D

    double maxd=0;
    for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) maxd=std::max(maxd, std::abs(VA(i,j)-VB(i,j)));
    EXPECT_NEAR(EhA, EhB, 1e-9);
    EXPECT_LT(maxd, 1e-9);
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
    cell.AddAtom(14, {0,0,0});                   // itsZ = TRUE species Z = 14 (Zion=4 supplied by the PP's ZionFn)
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D  lat(cell, ivec3_t(2,2,2));       // 2x2x2 = 8 k-points

    namespace L3=BasisSet::Lattice_3D;
    std::unique_ptr<BasisSet::Complex_BS> bs(L3::Factory(L3::Type::PW, lat, 4.0));
    BasisSet::irrepv_t irreps=bs->GetIrreps(Spin::None);   // one Bloch irrep per k-block (8)

    const int Nelec=8;
    Crystal_EC ec(irreps, Nelec);                          // Nval per k-block; weights handle the BZ sum

    cHamiltonian* ham=new Ham_PW_DFT(lat.GetStructure(), bs.get(), "Si", "LDA", 4);   // one-call: looks up + owns the GTH PP
    using qchem::SCFAccelerators::DIISParams;
    auto* acc=new qchem::SCFAccelerators::cSCFAcceleratorDIIS(DIISParams{8, 0.5, 1e-10, 1e-9});

    // Uniform-density seed on the first block: D=(N/n0)I gives rho(r)=N/V (uniform), the total density
    // every block's first Hartree/XC needs (a single block suffices since rho is constant).  Built
    // centrally by MakeSeedDensity (Uniform); also the plane-wave Default.
    qchem::SCFIterator::cSCFIterator scf(bs.get(), &ec, ham, acc,
                                         qchem::ChargeDensity::SeedStrategy::Uniform);

    SCFParams par;
    par.NMaxIter      =80;
    par.MinΔρ         =1e-7;   // tight: the FFT XC/Hartree made the multi-k run cheap, so complex DIIS now
                               // converges the BZ-sampled density fully (Etot=0.934176 exactly) in seconds.
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
    EXPECT_NEAR(charge,                  8.0,    1e-6);    // 8 valence electrons (BZ-weighted sum)
    EXPECT_NEAR(E.GetElectronicEnergy(), 0.934,  5e-3);    // band energy matches prototype ScfSiliconBZSampled
    // Physical total = electronic + ion-ion Ewald (same per-cell Madelung as Gamma) + G=0 alignment.
    EXPECT_NEAR(E.GetTotalEnergy(),     -7.7613, 5e-3) << "Enn="<<E.Enn<<" Ealign="<<E.Ealign;
    EXPECT_GT(gap, 0.0);                              // Si is a semiconductor
}

// Regression guard for the Hamiltonian-framework cache bug: rDynamic_HT_Imp::GetMatrix must invalidate
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

    std::unique_ptr<qchem::Hamiltonian::PW_Hartree> hart(NewPWHartree(F.pw));
    qchem::Hamiltonian::cDynamic_HT* ht=hart.get();
    ht->GetMatrix(&F.pw, Spin::None, &cd1);                       // populates the cache for cd1
    const chmat_t& M2 = ht->GetMatrix(&F.pw, Spin::None, &cd2);   // BUG: returns cd1's stale matrix

    double  Eh; chmat_t ref2=F.pw.Repulsion(cd2,Eh);             // the correct V_H for cd2
    double  diff=0.0;
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
            diff=std::max(diff, std::abs(dcmplx(M2(i,j))-dcmplx(ref2(i,j))));
    EXPECT_NEAR(diff, 0.0, 1e-10);   // fails today (M2 is cd1's matrix); passes after the cache-invalidation fix
}

} //namespace

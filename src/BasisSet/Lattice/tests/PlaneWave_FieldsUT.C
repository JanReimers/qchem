// File: src/BasisSet/Lattice/tests/PlaneWave_FieldsUT.C  Plane-wave basis / FFT / Poisson unit tests (moved out of
// IntegrationTests/PlaneWaveDFTUT.C on 2026-09-15, doc/TestSuitePlan.md phase 5): no Hamiltonian, no SCF.
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

import qchem.Tests.PW_Fields;   // the direct-grid oracles + PWFixture
import qchem.BasisSet.PlaneWave.PlaneWave_IBS;
import qchem.BasisSet.PlaneWave.Evaluators;   // PW_Grid_Evaluator -- a UNIT TEST may reach the internal
import qchem.BasisSet.Lattice.BasisSet;   // Factory(Type::PW, lat, Ecut, loc, nl) -> Complex_BS*
import qchem.ScalarFunction;                 // ScalarFunction<double> -- arg of the moved-here field oracles
import qchem.Mesh;                           // qcMesh::MeshParams (Vee_Hartree's fit-basis factory arg; ignored)
import qchem.Lattice_3D;     // UnitCell, Lattice_3D, ReciprocalLattice
import qchem.Ewald;          // EwaldEnergy (ion-ion Madelung term -> physical total energy)
import qchem.Types;          // dcmplx, ivec3_t, rvec_t, mat_t, chmat_t
import qchem.Blaze;          // mat_t<dcmplx>
import qchem.Math;           // Pi
import qchem.Pseudopotential.GTH_Potentials;    // GetGTH (CP2K GTH/HGH database reader)
import qchem.Fitting.FunctionFitter;                // Factory / ProjectedDensity_G / FunctionFitter_Density (item B)
import qchem.BasisSet.G_FieldEvaluator;             // the grid-engine seam (GridPoints/RhoOnGrid/Integral) for the item-K probe
import qchem.BasisSet.GMap;                             // ΔG_Map (the G-space coefficient map)
import qchem.BasisSet.Orbital_DFT_IBS;                      // cFIT_CD_ABS (the ortho density-fit basis face)
import qchem.Symmetry.Irrep;                        // Irrep
import qchem.LASolver;                              // complex Hermitian eigensolver
import qchem.Structure;                             // Molecule, Atom (the Si diamond basis)
import qchem.Matrix3D;                              // Matrix3D<double> (the FCC cell matrix)
import qchem.BasisSet.Internal.BasisSetImp;         // BasisSetImp<dcmplx> (single-block BasisSet container)

using namespace qchem;
using namespace qchem::tests::pw;

using BasisSet::PlaneWave::PlaneWave_IBS;
using BasisSet::PlaneWave::PW_Grid_Evaluator;   // internal grid evaluator (unit test may reach it directly)
using Pseudopotential::HGH_LocalPotential;
using Pseudopotential::HGH_SeparablePotential;
using Pseudopotential::GetGTH;
using Pseudopotential::GTH_PP;

namespace
{

// A single occupied plane wave (f=2 on band 0 = the i=0 PW) is a uniform density: the ONLY nonzero
// Fourier component is rho~(0)=N/Omega=2/Omega; every dm!=0 is exactly zero.
TEST(PlaneWaveFields, SinglePlaneWaveIsUniform)
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
TEST(PlaneWaveFields, ChargeSumRule)
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
TEST(PlaneWaveFields, HermitianSymmetry)
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
TEST(PlaneWaveFields, HartreeSingleCosineMatchesPoisson)
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


// The Hartree matrix from OverlapMatrix picks up V_H~(m_i-m_j) on each pair, and is Hermitian.
TEST(PlaneWaveFields, HartreeMatrixElementsAndHermiticity)
{
    PWFixture F;
    const double A=0.01;
    RhoG rho;
    rho[ivec3_t( 1,0,0)] = A;
    rho[ivec3_t(-1,0,0)] = A;
    double g02=(2*Pi/F.a)*(2*Pi/F.a);

    chmat_t V=F.pw.OverlapMatrix(HartreeVtilde(rho, F.B()));
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
TEST(PlaneWaveFields, DensityGridTransformRoundtrip)
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


// --- Stage 2: the basis-level high-level DFT capability, validated against the
// prototype's analytic/free-function results.  These are the questions the framework terms will ask --
// the term hands a real-space ScalarFunction and the basis owns the integration (no G-vectors exposed).

// Integral(f) = integral f d3r over the cell: a constant integrates to const*Omega; a reciprocal
// cosine integrates to zero.
TEST(PlaneWaveFields, BasisIntegralScalar)
{
    PWFixture F;
    FieldFn cst([](const rvec3_t&){return 2.5;});
    EXPECT_NEAR(IntegralField(F.pw, cst), 2.5*F.Omega, 1e-9*F.Omega);

    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    FieldFn cosfn([&](const rvec3_t& r){return std::cos(G0*r);});
    EXPECT_NEAR(IntegralField(F.pw, cosfn), 0.0, 1e-9*F.Omega);
}


// Repulsion on a single-cosine density rho = rho0 + 2A cos(G0.r) reproduces the Poisson result
// (E_H = Omega 4 pi A^2/|G0|^2, V_H~(G0) = 4 pi A/|G0|^2) -- same check as HartreeSingleCosineMatchesPoisson,
// now driven entirely through the basis's high-level method (the basis samples + transforms internally).
TEST(PlaneWaveFields, BasisRepulsionMatchesPoisson)
{
    PWFixture F;
    const double rho0=2.0/F.Omega, A=0.01;
    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    double  g02=G0*G0;
    FieldFn rho([&](const rvec3_t& r){return rho0 + 2*A*std::cos(G0*r);});

    chmat_t VH=RepulsionField(F.pw, rho);   // Hartree matrix; the E_H-vs-Poisson check is HartreeSingleCosineMatchesPoisson

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
TEST(PlaneWaveFields, BasisIntegralPotentialConstAndCosine)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();

    FieldFn cst([](const rvec3_t&){return 0.7;});
    chmat_t Vc=OverlapField(F.pw, cst);
    for (size_t i=0;i<n;i++)
    {
        EXPECT_NEAR(std::real(dcmplx(Vc(i,i))), 0.7, 1e-9);
        for (size_t j=i+1;j<n;j++) EXPECT_NEAR(std::abs(dcmplx(Vc(i,j))), 0.0, 1e-9);
    }

    rvec3_t G0=F.B().ToCartesian(rvec3_t(1,0,0));
    FieldFn cosfn([&](const rvec3_t& r){return std::cos(G0*r);});
    chmat_t Vk=OverlapField(F.pw, cosfn);
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


// Item B: the ortho (plane-wave) density fitter's fitted field is a REAL, evaluatable ScalarFunction --
// rho_fit(r) = Re Σ_dm rho~(dm) e^{i(B·dm)·r} (what the GUI plots).  Fit a non-uniform density through the
// production Factory path, then check the fitter's op(r)/Gradient via the G_FieldEvaluator DIP seam against
// (a) an independent inverse transform, (b) the FFT route RhoOnGrid at r=0, (c) finite-difference gradient.
TEST(PlaneWaveFields, OrthoFitterRealSpaceField)
{
    PWFixture F;
    size_t n=F.pw.GetNumFunctions();
    hmat_t<dcmplx> D=blazem::zeroH<dcmplx>(n);
    D(0,0)=1.0; D(0,1)=1.0; D(1,1)=1.0;                   // Hermitian, non-uniform -> a non-trivial rho~(dm)
    ΔG_Map rhoTilde=RhoTilde(F.pw, D);

    // The Factory-built ortho density fitter (the production path).
    auto fb=std::shared_ptr<const qchem::BasisSet::cFIT_CD_ABS>(F.pw.CreateCDFitBasisSet(nullptr, qcMesh::MeshParams{}));
    auto fitter=qchem::Fitting::Factory(fb);
    fitter->DoFit(qchem::Fitting::ProjectedDensity_G(rhoTilde));
    // The fitted FIELD is a CAPABILITY, not part of the fitter face (2026-08-24): an AO fit could only
    // deliver it by handing its coefficients back to the basis, so the whole chain went.  A plane-wave fit
    // genuinely can -- it inverse-transforms its own {G} -- so it derives ScalarFunction itself and a
    // consumer asks.  Same cross-cast the scalar-fitter sibling above already makes.
    const ScalarFunction<double>& rhoFit=dynamic_cast<const ScalarFunction<double>&>(*fitter);

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
    rvec_t grid=GridOf(F.pw).RhoOnGrid(rhoTilde);
    EXPECT_NEAR(rhoFit(rvec3_t(0.0,0.0,0.0)), grid[0], 1e-11);

    // Gradient: finite-difference the field at a generic point.
    rvec3_t r0(1.3,0.7,2.1), g=rhoFit.Gradient(r0);
    const double h=1e-5;
    EXPECT_NEAR(g.x, (rhoFit(rvec3_t(r0.x+h,r0.y,r0.z))-rhoFit(rvec3_t(r0.x-h,r0.y,r0.z)))/(2*h), 1e-4);
    EXPECT_NEAR(g.y, (rhoFit(rvec3_t(r0.x,r0.y+h,r0.z))-rhoFit(rvec3_t(r0.x,r0.y-h,r0.z)))/(2*h), 1e-4);
    EXPECT_NEAR(g.z, (rhoFit(rvec3_t(r0.x,r0.y,r0.z+h))-rhoFit(rvec3_t(r0.x,r0.y,r0.z-h)))/(2*h), 1e-4);
}

} //namespace

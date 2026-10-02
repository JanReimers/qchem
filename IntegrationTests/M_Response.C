// File: IntegrationTests/M_Response.C  Molecular linear response through the qchem::Calculation facade
// (doc/LinearResponsePlan.md stage R1: CPHF static polarisability).
//
// THE ORACLE is PySCF on the SAME molecule, geometry and basis (scripts/r1_h2o_polarizability.py, banked
// 2026-09-28): dzvp Cartesian, whose HF energy matches ours to 1e-9 (M_Calculation.WaterEnergy), so the HF
// polarisability is a TIGHT gate -- analytic CPHF and finite field agree there.
#include "gtest/gtest.h"
#include <cmath>
#include <stdexcept>

import qchem.Calculation;
import qchem.StructureData;        // GetMolecule("H2O") -- the shared water geometry (D-MAKEWATER)
import qchem.Structure;
import qchem.Types;        // Vector3D
import qchem.Blaze;
import qchem.SCFParams;
using namespace qchem;


// PySCF RHF/dzvp(cart) CPHF, scripts/r1_h2o_polarizability.py: diag(3.19770284, 7.11975051, 5.54545175) bohr^3.
static const double alphaPySCF[3]={3.19770284, 7.11975051, 5.54545175};

// A response linearises about a CONVERGED state (§4c): tight enough that the gates below measure the response.
static const SCFParams tight = {.NMaxIter=80, .MinΔρ=1e-9, .MinΔFD=1e-10, .MinVirial=1e2};

static rmat_t Alpha(const Calculation& calc)
{
    auto a=calc.StaticPolarizability();
    EXPECT_TRUE(a.IsOk()) << (a ? "" : a.Error().detail);
    return a ? a.Value() : rmat_t(3,3,0.0);
}

TEST(M_Response, HF_Water_PolarizabilityIsPySCF)
{
    Calculation calc(qchem::StructureData::GetMolecule("H2O"), {.basis="dzvp"});
    ASSERT_TRUE(calc.Converge(tight));
    const rmat_t a=Alpha(calc);
    for (size_t i=0;i<3;i++)
    {
        EXPECT_NEAR(a(i,i), alphaPySCF[i], 1e-6*alphaPySCF[i]) << "alpha_" << i << i;
        for (size_t j=0;j<3;j++) if (j!=i) EXPECT_NEAR(a(i,j), 0.0, 1e-6) << "alpha_" << i << j << " (C2v: diagonal)";
    }
}

//! Pol is the PRIMARY formulation (CLAUDE.md): the same closed shell imposed POLARIZED -- per-channel exchange,
//! spin-resolved transition density -- must give the same polarisability as the folded doublet.
TEST(M_Response, HF_Water_Pol_eqUnPol)
{
    Calculation unpol(qchem::StructureData::GetMolecule("H2O"), {.basis="dzvp"});
    Calculation pol  (qchem::StructureData::GetMolecule("H2O"), {.basis="dzvp", .spin=SpinGroup::Polarized});
    ASSERT_TRUE(unpol.Converge(tight));
    ASSERT_TRUE(pol  .Converge(tight));
    const rmat_t a=Alpha(unpol), b=Alpha(pol);
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) EXPECT_NEAR(a(i,j), b(i,j), 1e-8) << "alpha_" << i << j;
}

//! Ruling Q5: the dipole matrices are NUMERICAL, so their error is MEASURED -- a finer mesh must not move alpha
//! beyond the floor the Becke quadrature sets.  MEASURED 2026-09-28 (bohr^3, xx / yy / zz):
//!     default  MHL 80 / Lebedev 35      3.1977039  7.1197506  5.5454507
//!              MHL 120 / GaussLeg 47    3.1977031  7.1197503  5.5454510
//!              MHL 250 / GaussLeg 71    3.1977026  7.1197503  5.5454511
//!     PySCF (analytic <r>)              3.1977028  7.1197505  5.5454518
//! i.e. the numerical dipole OSCILLATES about the analytic answer at ~1e-7 relative: a quadrature floor, not a
//! trend.  Gated at 1e-6 relative (4x the largest move seen).  ⚠ If a gate ever needs better than ~1e-7, the
//! cure is ANALYTIC dipole integrals (user, 2026-09-28: MnD Hermite set-up, libcint int1e_r as the oracle), not
//! a bigger mesh.
TEST(M_Response, HF_Water_DipoleMeshConverged)
{
    Calculation calc(qchem::StructureData::GetMolecule("H2O"), {.basis="dzvp"});
    ASSERT_TRUE(calc.Converge(tight));
    const rmat_t a=Alpha(calc);
    auto fine=calc.StaticPolarizability({.radial=qcMesh::RadialKind::MHL, .nRadial=120, .mhl_m=3, .mhl_alpha=2.0,
                                         .angular=qcMesh::AngularKind::GaussLegendre, .angularDegree=47, .beckeOrder=3});
    ASSERT_TRUE(fine.IsOk()) << fine.Error().detail;
    for (size_t i=0;i<3;i++) EXPECT_NEAR(a(i,i), fine.Value()(i,i), 1e-6*a(i,i)) << "alpha_" << i << i << " moved with the dipole mesh";
}

//! A SYMMETRY-ADAPTED molecule is REFUSED, loudly: x and y are not totally symmetric in C2v, so they couple
//! different irreps, and without the point-group product selection rule their response would be silently zero.
TEST(M_Response, HF_Water_SymmetryAdaptedIsRefused)
{
    Calculation calc(qchem::StructureData::GetMolecule("H2O"), {.basis="dzvp", .symmetry=true});
    EXPECT_THROW((void)calc.StaticPolarizability(), std::logic_error);
}

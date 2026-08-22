// File: RealComplexTermsUT.C  Step-3c-1 gates: the periodic terms serve a REAL TRIM block (doc/RealComplexPlan.md).
//
// Step 3a made the real TRIM block (tGPW_IBS<double>) exist; this increment makes the TERM STACK able to
// answer it.  Each periodic term now carries the Static/Dynamic_HT_RealBlock capability face (the
// V1.6/V1.7 cross-cast idiom), implemented by ONE scalar-generic assembly body per term -- so a real
// block's matrix is the SAME arithmetic narrowed at the end, and must equal the complex block's real
// part BITWISE (EXPECT_EQ on doubles; Step 0's exact +/-1 phases are what make zero tolerance possible).
//
// Covered here: the static set (Kinetic, Ven_PP_Short/Long/NonLocal) and the density-dependent
// Vee_Hartree + the Becke-route Vxc_Quadrature.  Each scalar arm drives its OWN production-shaped
// density (a real child on the real block, a complex child on the complex block -- the 3c-3 mixed
// shape): the GPW integrate-back's D-aware screen is INSTANCE-PAIRED state (the collocation and its
// adjoint must go through the same block instance -- the 2026-08-18 cross-run-pollution fix scoped the
// 3C tensors per instance, see tGPW_IBS::Repulsion3C), so a cross-asked arrangement (density collocated
// on cx, matrix asked of re) is no longer screen-consistent and was only ever bitwise by the two twins
// accidentally sharing one DBCache'd closure.  The pair-quadrature RAW-route real path needs the full
// collocation feed and is gated with the 3c-2 SCF wiring.
#include "gtest/gtest.h"
#include <complex>
#include <memory>

import qchem.Structure;
import qchem.UnitCell;
import qchem.BasisSet;
import qchem.BasisSet.Molecule.Factory;
import qchem.BasisSet.Lattice_3D.GPW_IBS;          // tGPW_IBS<T> (both block alternatives)
import qchem.Hamiltonian;                          // Static/Dynamic_HT_RealBlock (the capability faces)
import qchem.Hamiltonian.Internal.Kinetic;         // Kinetic<T> (tests may cheat-import internals)
import qchem.Hamiltonian.Internal.Hamiltonian;     // tHamiltonianImp<dcmplx> (the 3c-2 assembly gate)
import qchem.Hamiltonian.Internal.PWTerms;         // the periodic term set + the two XC_Quadrature strategies
import qchem.BasisSet.DeltaFit_IBS;                // DeltaFit_IBS -- the delta basis the singles strategy runs on
import qchem.Symmetry.Factory;                     // BlochFactory (its Gamma irrep)
import qchem.Hamiltonian.Internal.SlaterExchange;  // SlaterExchange (the Dirac-exchange functional)
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH (Si GTH-LDA-q4)
import qchem.ChargeDensity.Imp.IrrepCD;            // PeriodicIrrepCD (the periodic leaf -- 3c-2b)
import qchem.CompositeCD;                          // tComposite_CD<dcmplx> (the run-shaped density)
import qchem.Mesh;                                 // qcMesh::MeshParams
import qchem.Fitting.FunctionFitter;               // CoulombMetric_ProjectedDensity (the honest-probe check)
import qchem.Blaze;
import qchem.Types;

using namespace qchem;
using namespace qchem::Hamiltonian;
using BasisSet::Real_BS;
using BasisSet::Lattice_3D::tGPW_IBS;

namespace
{
struct Rig
{
    UnitCell cell;
    std::shared_ptr<const Structure> st;
    std::shared_ptr<const Real_BS>   mol;
    std::unique_ptr<tGPW_IBS<double>> re;   // the real Gamma TRIM block (Step 3a)
    std::unique_ptr<tGPW_IBS<dcmplx>> cx;   // the same physical block, complex
    Rig() : cell(8.0)
    {
        cell.AddAtom(14,{0.5,0.5,0.5});
        st=cell.Clone();
        namespace M = qchem::BasisSet::Molecule;
        mol.reset(M::Factory(M::BasisSetData::SIPP, &cell, M::Engine::MnD, M::Angular::Cartesian));
        re=std::make_unique<tGPW_IBS<double>>(cell, ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
        cx=std::make_unique<tGPW_IBS<dcmplx>>(cell, ivec3_t(2,2,2), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    }
    //! A one-block complex density in the run's shape (tComposite_CD<dcmplx>): D(0,0)=2 -- two electrons
    //! in the first orbital.  Enough to make V_H / rho genuinely nonzero; the terms never see past the
    //! composite face.
    std::unique_ptr<ChargeDensity::cDM_CD> MakeDensity() const
    {
        using namespace qchem::ChargeDensity;
        chmat_t D(cx->GetNumFunctions());
        D(0,0)=2.0;
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(
            new PeriodicIrrepCD<dcmplx>(D, cx.get(), cx->GetIrrep(Spin::None))));
        return cd;
    }
    //! The REAL twin of \c MakeDensity: the same D(0,0)=2 as a REAL child on the REAL block inside the
    //! complex-faced composite -- the 3c-3 mixed shape.  The real arm of every density-dependent gate
    //! drives THIS density, so its collocation and integrate-back pair on the real instance exactly as
    //! a production mixed run pairs them (the screen-consistency note in the file header).
    std::unique_ptr<ChargeDensity::cDM_CD> MakeRealDensity() const
    {
        using namespace qchem::ChargeDensity;
        rsmat_t D(re->GetNumFunctions());
        D(0,0)=2.0;
        auto cd=std::make_unique<tComposite_CD<dcmplx>>();
        cd->Insert(std::unique_ptr<tDM_CD<double>>(
            new PeriodicIrrepCD<double>(D, re.get(), re->GetIrrep(Spin::None))));
        return cd;
    }
};

void ExpectBitwiseEqual(const rsmat_t& R, const chmat_t& C, const char* what)
{
    ASSERT_EQ(R.rows(), C.rows()) << what;
    for (size_t i=0;i<R.rows();i++)
        for (size_t j=0;j<R.columns();j++)
            EXPECT_EQ(R(i,j), std::real(C(i,j))) << what << " (" << i << "," << j << ")";
}
// The Becke quadrature GEMM comparison CANNOT be bitwise, and that is the Step-3 win showing itself:
// the real block's Phi table drives blaze's REAL dgemm-shaped kernel, whose accumulation order differs
// from the complex kernel's -- same exact integrand, last-ulp summation differences (~3e-15 measured).
// Everything ELEMENTWISE or same-order (S/T/V/KB, the Hartree ContractAdjoint) stays EXPECT_EQ above.
void ExpectMachineEqual(const rsmat_t& R, const chmat_t& C, const char* what)
{
    ASSERT_EQ(R.rows(), C.rows()) << what;
    for (size_t i=0;i<R.rows();i++)
        for (size_t j=0;j<R.columns();j++)
            EXPECT_NEAR(R(i,j), std::real(C(i,j)), 1e-13) << what << " (" << i << "," << j << ")";
}
} //anon

// The STATIC set: kinetic + the three PP terms.  One capability cross-cast each, then bitwise equality.
TEST(RealComplexTerms, StaticTermsServeTheRealBlockBitwise)
{
    Rig rig;
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    Kinetic<dcmplx> kin;
    Ven_PP_Short    shrt(rig.st, &gth.local);
    Ven_PP_Long     lng (rig.st, &gth.local);
    Ven_PP_NonLocal nl  (rig.st, &gth.nonlocal);

    const std::vector<std::pair<const cStatic_HT*,const char*>> terms =
        {{&kin,"Kinetic"},{&shrt,"Ven_PP_Short"},{&lng,"Ven_PP_Long"},{&nl,"Ven_PP_NonLocal"}};
    for (auto [t,name] : terms)
    {
        auto* rb=dynamic_cast<const Static_HT_RealBlock*>(t);
        ASSERT_NE(rb,nullptr) << name << " must carry the real-block capability (Step 3c)";
        ExpectBitwiseEqual(rb->GetMatrix(rig.re.get(), Spin::None),
                           t ->GetMatrix(rig.cx.get(), Spin::None), name);
    }
}

// The DYNAMIC pair reachable without the full SCF: Hartree (G-space Poisson off the shared V_H(G)) and
// the Becke-route Vxc_Quadrature (the engine's rho raster is block-independent; the real block gets its
// own REAL Phi table, so its quadrature GEMM runs in real arithmetic -- the first realized Step-3 win).
TEST(RealComplexTerms, HartreeAndBeckeXcServeTheRealBlockBitwise)
{
    Rig rig;
    auto cd  = rig.MakeDensity();       // the complex arm's density (complex child on cx)
    auto cdr = rig.MakeRealDensity();   // the real arm's TWIN (real child on re -- the 3c-3 mixed shape)

    Vee_Hartree::fbs_t fb(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    Vee_Hartree vh(fb);
    {
        auto* rb=dynamic_cast<const Dynamic_HT_RealBlock*>(&vh);
        ASSERT_NE(rb,nullptr) << "Vee_Hartree must carry the real-block capability (Step 3c)";
        // Concrete-typed calls are ambiguous across the two GetMatrix base classes; go through the faces.
        const chmat_t Vc=static_cast<const cDynamic_HT&>(vh).GetMatrix(rig.cx.get(), Spin::None, cd.get());   // complex first
        const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, cdr.get());
        ExpectBitwiseEqual(Vr, Vc, "Vee_Hartree");
    }

    auto q = rig.cx->CreateXCQuadrature(rig.st.get(), qcMesh::MeshParams{});
    auto dfb=std::make_shared<const BasisSet::DeltaFit_IBS>(q, Symmetry::BlochFactory(ivec3_t(1,1,1), ivec3_t(0,0,0)));
    auto engine=std::make_shared<const XC_SinglesQuadrature>(dfb);
    Vxc_Quadrature vxc(std::make_shared<SlaterExchange>(2.0/3.0), engine);
    {
        auto* rb=dynamic_cast<const Dynamic_HT_RealBlock*>(&vxc);
        ASSERT_NE(rb,nullptr) << "Vxc_Quadrature must carry the real-block capability (Step 3c)";
        const chmat_t Vc=static_cast<const cDynamic_HT&>(vxc).GetMatrix(rig.cx.get(), Spin::None, cd.get());
        const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, cdr.get());
        ExpectMachineEqual(Vr, Vc, "Vxc_Quadrature");
    }
}


// PRODUCTION-SHAPED rig (Step 3c-3 divergence hunt): the FCC-diamond 2-atom cell the acceptance runs,
// where the single-atom cubic Rig above cannot reach multi-centre / oblique-cell defects.
struct FccRig
{
    FCCUnitCell cell;
    std::shared_ptr<const Structure> st;
    std::shared_ptr<const Real_BS>   mol;
    std::unique_ptr<tGPW_IBS<double>> re;
    std::unique_ptr<tGPW_IBS<dcmplx>> cx;
    FccRig() : cell(10.26)
    {
        cell.AddAtom(14,{0,0,0});
        cell.AddAtom(14,{0.25,0.25,0.25});
        st=cell.Clone();
        namespace M = qchem::BasisSet::Molecule;
        mol.reset(M::Factory(M::BasisSetData::SIPP_SR, &cell, M::Engine::MnD, M::Angular::Cartesian));
        re=std::make_unique<tGPW_IBS<double>>(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
        cx=std::make_unique<tGPW_IBS<dcmplx>>(cell, ivec3_t(1,1,1), ivec3_t(0,0,0), mol, /*densityEcut*/20.0);
    }
};

TEST(RealComplexTerms, StaticTermsServeTheRealBlockBitwise_FccDiamond)
{
    FccRig rig;
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    Kinetic<dcmplx> kin;
    Ven_PP_Short    shrt(rig.st, &gth.local);
    Ven_PP_Long     lng (rig.st, &gth.local);
    Ven_PP_NonLocal nl  (rig.st, &gth.nonlocal);

    const std::vector<std::pair<const cStatic_HT*,const char*>> terms =
        {{&kin,"Kinetic"},{&shrt,"Ven_PP_Short"},{&lng,"Ven_PP_Long"},{&nl,"Ven_PP_NonLocal"}};
    for (auto [t,name] : terms)
    {
        auto* rb=dynamic_cast<const Static_HT_RealBlock*>(t);
        ASSERT_NE(rb,nullptr) << name;
        ExpectBitwiseEqual(rb->GetMatrix(rig.re.get(), Spin::None),
                           t ->GetMatrix(rig.cx.get(), Spin::None), name);
    }
}

// The SEED-SHAPED density on the production cell (3c-3 hunt): the run's two arms feed the SAME uniform
// D = (N/n)·I through DIFFERENT leaves (real vs complex).  Their collocated rho and V_H(G) must agree
// bitwise, or the very first Fock differs.
TEST(RealComplexTerms, SeedDensityTrioMatches_FccDiamond)
{
    using namespace qchem::ChargeDensity;
    FccRig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); chmat_t Dc(n);
    for (size_t i=0;i<n;i++) { Dr(i,i)=8.0/double(n); Dc(i,i)=8.0/double(n); }
    PeriodicIrrepCD<double> cdr(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None));
    PeriodicIrrepCD<dcmplx> cdc(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None));
    EXPECT_NEAR(cdr.GetTotalCharge(), cdc.GetTotalCharge(), 1e-13);

    Vee_Hartree::fbs_t fb(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> sf(rig.cx->CreateVxcFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    ΔG_Map rr=cdr.GetRepulsion3C(*fb), rc=cdc.GetRepulsion3C(*fb);
    ASSERT_EQ(rr.size(), rc.size());
    for (const auto& [dm,v] : rr)
    {
        auto it=rc.find(dm);
        ASSERT_NE(it, rc.end());
        EXPECT_EQ(v, it->second) << "V_H(dm) mismatch";
    }
    const rvec_t gr=cdr.GetRhoOnGrid(*sf), gc=cdc.GetRhoOnGrid(*sf);
    ASSERT_EQ(gr.size(), gc.size());
    for (size_t g=0; g<gr.size(); g++) EXPECT_EQ(gr[g], gc[g]) << "raw rho_DM mismatch at g="<<g;

    // ...and the SAME through the run-shaped COMPOSITES (the exact seed objects the two arms feed the
    // raw XC route): a mixed composite's rho visit must reach the real child identically.
    tComposite_CD<dcmplx> mixed;  mixed.Insert(std::unique_ptr<tDM_CD<double>>(
        new PeriodicIrrepCD<double>(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None))));
    tComposite_CD<dcmplx> cplx;   cplx.Insert(std::unique_ptr<tDM_CD<dcmplx>>(
        new PeriodicIrrepCD<dcmplx>(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None))));
    const auto& fm=dynamic_cast<const FourierDensity&>(mixed);
    const auto& fx=dynamic_cast<const FourierDensity&>(cplx);
    const rvec_t cgr=fm.GetRhoOnGrid(*sf), cgc=fx.GetRhoOnGrid(*sf);
    ASSERT_EQ(cgr.size(), cgc.size());
    for (size_t g=0; g<cgr.size(); g++) EXPECT_EQ(cgr[g], cgc[g]) << "composite raw rho mismatch at g="<<g;
}

// The RAW-route pair-quadrature gate (Step 3c-3): the production Uniform-mesh XC.  The collocated rho_DM
// feed is block-independent; the real block's matrix comes back through the SAME raw adjoint tensor
// (TFit==dcmplx either way) narrowed at the end, so it must equal the complex block's real part
// BITWISE.  (This was the one term route 3c-1 could not gate without the full collocation feed.)
TEST(RealComplexTerms, RawRouteXcServesTheRealBlockBitwise)
{
    Rig rig;
    auto cd  = rig.MakeDensity();       // the complex arm's density
    auto cdr = rig.MakeRealDensity();   // the real arm's twin (screen-consistent pairing; file header)
    XC_PairQuadrature::fbs_t fb(rig.cx->CreateVxcFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    Vxc_Quadrature vxc(std::make_shared<SlaterExchange>(2.0/3.0), MakeXCQuadrature(fb));
    auto* rb=dynamic_cast<const Dynamic_HT_RealBlock*>(&vxc);
    ASSERT_NE(rb,nullptr) << "Vxc_Quadrature must carry the real-block capability (Step 3c)";
    const chmat_t Vc=static_cast<const cDynamic_HT&>(vxc).GetMatrix(rig.cx.get(), Spin::None, cd.get());
    const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, cdr.get());
    ExpectBitwiseEqual(Vr, Vc, "Vxc_Quadrature (pair) raw");
}

// The ASSEMBLY gate (Step 3c-2): the Hamiltonian's Ham_RealBlock fold over the term faces must equal the
// native complex assembly's real part BITWISE -- each term's block already does (the 3c-1 gates), and the
// fold accumulates them elementwise in the same order.  This is the exact matrix a tIrrepWF<double> child
// receives from its CalculateH inside a complex run.
namespace { struct MixedHam : Hamiltonian::tHamiltonianImp<dcmplx> {}; }

TEST(RealComplexTerms, HamiltonianAssemblyServesTheRealBlockBitwise)
{
    Rig rig;
    auto cd  = rig.MakeDensity();       // the complex arm's density
    auto cdr = rig.MakeRealDensity();   // the real arm's twin (screen-consistent pairing; file header)
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    MixedHam ham;
    ham.Add(new Kinetic<dcmplx>);
    ham.Add(new Ven_PP_Short   (rig.st, &gth.local));
    ham.Add(new Ven_PP_Long    (rig.st, &gth.local));
    ham.Add(new Ven_PP_NonLocal(rig.st, &gth.nonlocal));
    ham.Add(new Vee_Hartree(Vee_Hartree::fbs_t(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}))));

    const chmat_t Hc = ham.GetMatrix(rig.cx.get(), Spin::None, cd.get(), nullptr);
    auto& rb = dynamic_cast<Hamiltonian::Ham_RealBlock&>(ham);   // the face a real WF child drives
    const rsmat_t Hr = rb.GetMatrix(rig.re.get(), Spin::None, cdr.get(), nullptr);
    ExpectBitwiseEqual(Hr, Hc, "Ham_RealBlock fold");
}


// ===== Step 3c-2b: the PERIODIC LEAF SPLIT + the mixed-composite energy/density arms ==================
// PeriodicIrrepCD<double> is the real TRIM block's density: real D, real DM-side arithmetic, and the
// SAME reciprocal trio as its complex twin -- the lineage-as-class split (V1.7's scalar⇔lineage
// conflation retired).  The trio must match the complex twin BITWISE: the contraction runs the same
// tensor in the same order, only D's (exactly zero) imaginary parts differ.
TEST(RealComplexTerms, RealPeriodicLeafFourierTrioMatchesComplexBitwise)
{
    using namespace qchem::ChargeDensity;
    Rig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); Dr(0,0)=2.0;
    chmat_t Dc(n); Dc(0,0)=2.0;
    PeriodicIrrepCD<double> cdr(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None));
    PeriodicIrrepCD<dcmplx> cdc(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None));

    // Both leaves carry the Fourier face; NEITHER carries the finite (AO/HF) faces -- the honest probes.
    EXPECT_NE(dynamic_cast<const FourierDensity*>(&cdr), nullptr);
    EXPECT_EQ(dynamic_cast<const Fitting::CoulombMetric_ProjectedDensity*>(&cdr), nullptr);
    EXPECT_EQ(dynamic_cast<const tHF_Pair_CD<double>*>(&cdr), nullptr);

    Vee_Hartree::fbs_t fb(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    std::shared_ptr<const BasisSet::cFIT_SF_ABS> sf(rig.cx->CreateVxcFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    ΔG_Map rr=cdr.GetRepulsion3C(*fb), rc=cdc.GetRepulsion3C(*fb);
    ASSERT_EQ(rr.size(), rc.size());
    for (const auto& [dm,v] : rr)
    {
        auto it=rc.find(dm);
        ASSERT_NE(it, rc.end());
        EXPECT_EQ(v, it->second) << "V_H(dm) mismatch";
    }
    const rvec_t gr=cdr.GetRhoOnGrid(*sf), gc=cdc.GetRhoOnGrid(*sf);
    ASSERT_EQ(gr.size(), gc.size());
    for (size_t g=0; g<gr.size(); g++) EXPECT_EQ(gr[g], gc[g]) << "raw rho_DM mismatch at g="<<g;
}

// The MIXED composite's energy/density arms: ONE real child inside a complex-faced composite must give
// the same static/dynamic energy contractions and the same rho(points) as the all-complex twin.
TEST(RealComplexTerms, MixedCompositeEnergyAndRhoMatchComplex)
{
    using namespace qchem::ChargeDensity;
    Rig rig;
    const size_t n=rig.re->GetNumFunctions();
    rsmat_t Dr(n); Dr(0,0)=2.0;
    chmat_t Dc(n); Dc(0,0)=2.0;

    tComposite_CD<dcmplx> mixed;   // complex face, REAL child (the 3c-3 shape)
    mixed.Insert(std::unique_ptr<tDM_CD<double>>(
        new PeriodicIrrepCD<double>(Dr, rig.re.get(), rig.re->GetIrrep(Spin::None))));
    tComposite_CD<dcmplx> cplx;    // the all-complex twin
    cplx.Insert(std::unique_ptr<tDM_CD<dcmplx>>(
        new PeriodicIrrepCD<dcmplx>(Dc, rig.cx.get(), rig.cx->GetIrrep(Spin::None))));

    // STATIC arm: the term's real-block face IS a tStatic_CC<double>, so the real child contracts natively.
    Kinetic<dcmplx> kin;
    EXPECT_EQ(mixed.DM_Contract(static_cast<const cStatic_CC*>(&kin)),
              cplx .DM_Contract(static_cast<const cStatic_CC*>(&kin)));

    // DYNAMIC arm: the run-typed GetEMatrixR client (E=D.V identity over the real cache).
    auto cdrun = rig.MakeDensity();   // the run's density (defines V_H)
    Vee_Hartree vh(Vee_Hartree::fbs_t(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{})));
    EXPECT_EQ(mixed.DM_Contract(static_cast<const cDynamic_CC*>(&vh), cdrun.get()),
              cplx .DM_Contract(static_cast<const cDynamic_CC*>(&vh), cdrun.get()));

    // RHO arm: the cross-scalar child self-evaluates pointwise; same rho(r) either way.
    rvec3vec_t pts(3);
    pts[0]=rig.cell.ToCartesian(rvec3_t(0.3,0.4,0.7));
    pts[1]=rig.cell.ToCartesian(rvec3_t(0.5,0.5,0.5));
    pts[2]=rig.cell.ToCartesian(rvec3_t(0.1,0.9,0.2));
    const rvec_t rm=mixed.DM_RhoAtPoints(pts, {});
    const rvec_t rc=cplx .DM_RhoAtPoints(pts, {});
    for (size_t g=0; g<pts.size(); g++) EXPECT_EQ(rm[g], rc[g]) << "rho mismatch at point "<<g;

    // And the composite Fourier trio (the generic visit) now reaches the REAL child's face too.
    // Through the FourierDensity FACE: the composite's own AO-projection GetRepulsion3C(rFIT*) hides the
    // Fourier overload at concrete type (member-name hiding); the face is how every consumer reaches it.
    Vee_Hartree::fbs_t fb2(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    const auto& fm=dynamic_cast<const FourierDensity&>(mixed);
    const auto& fx=dynamic_cast<const FourierDensity&>(cplx);
    ΔG_Map vm=fm.GetRepulsion3C(*fb2), vc=fx.GetRepulsion3C(*fb2);
    ASSERT_EQ(vm.size(), vc.size());
    for (const auto& [dm,v] : vm) { auto it=vc.find(dm); ASSERT_NE(it,vc.end()); EXPECT_EQ(v,it->second); }
}

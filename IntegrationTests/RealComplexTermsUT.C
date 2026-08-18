// File: RealComplexTermsUT.C  Step-3c-1 gates: the periodic terms serve a REAL TRIM block (doc/RealComplexPlan.md).
//
// Step 3a made the real TRIM block (tGPW_IBS<double>) exist; this increment makes the TERM STACK able to
// answer it.  Each periodic term now carries the Static/Dynamic_HT_RealBlock capability face (the
// V1.6/V1.7 cross-cast idiom), implemented by ONE scalar-generic assembly body per term -- so a real
// block's matrix is the SAME arithmetic narrowed at the end, and must equal the complex block's real
// part BITWISE (EXPECT_EQ on doubles; Step 0's exact +/-1 phases are what make zero tolerance possible).
//
// Covered here: the static set (Kinetic, Ven_PP_Short/Long/NonLocal) and the density-dependent
// Vee_Hartree + the Becke-route DeltaFittedVxc (driven by a hand-built one-block complex density -- the
// term's density-dependent state is block-independent, which is exactly what these gates pin).  The
// PWFittedVxc RAW-route real path needs the full collocation feed and is gated with the 3c-2 SCF wiring.
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
import qchem.Hamiltonian.Internal.PWTerms;         // the periodic term set + XC_GridEngine
import qchem.Hamiltonian.Internal.SlaterExchange;  // SlaterExchange (the Dirac-exchange functional)
import qchem.Pseudopotential.GTH_Potentials;       // GetGTH (Si GTH-LDA-q4)
import qchem.ChargeDensity.Imp.IrrepCD;            // IrrepCD<dcmplx> (hand-built block density)
import qchem.CompositeCD;                          // tComposite_CD<dcmplx> (the run-shaped density)
import qchem.Mesh;                                 // qcMesh::MeshParams
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
            new IrrepCD<dcmplx>(D, cx.get(), cx->GetIrrep(Spin::None))));
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
// the Becke-route DeltaFittedVxc (the engine's rho raster is block-independent; the real block gets its
// own REAL Phi table, so its quadrature GEMM runs in real arithmetic -- the first realized Step-3 win).
TEST(RealComplexTerms, HartreeAndBeckeXcServeTheRealBlockBitwise)
{
    Rig rig;
    auto cd = rig.MakeDensity();

    Vee_Hartree::fbs_t fb(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}));
    Vee_Hartree vh(fb);
    {
        auto* rb=dynamic_cast<const Dynamic_HT_RealBlock*>(&vh);
        ASSERT_NE(rb,nullptr) << "Vee_Hartree must carry the real-block capability (Step 3c)";
        // Concrete-typed calls are ambiguous across the two GetMatrix base classes; go through the faces.
        const chmat_t Vc=static_cast<const cDynamic_HT&>(vh).GetMatrix(rig.cx.get(), Spin::None, cd.get());   // complex first
        const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, cd.get());
        ExpectBitwiseEqual(Vr, Vc, "Vee_Hartree");
    }

    auto q = rig.cx->CreateXCQuadrature(rig.st.get(), qcMesh::MeshParams{});
    auto engine=std::make_shared<const XC_GridEngine>(q.mesh, q.fold);
    DeltaFittedVxc vxc(std::make_shared<SlaterExchange>(2.0/3.0), engine);
    {
        auto* rb=dynamic_cast<const Dynamic_HT_RealBlock*>(&vxc);
        ASSERT_NE(rb,nullptr) << "DeltaFittedVxc must carry the real-block capability (Step 3c)";
        const chmat_t Vc=static_cast<const cDynamic_HT&>(vxc).GetMatrix(rig.cx.get(), Spin::None, cd.get());
        const rsmat_t Vr=rb->GetMatrix(rig.re.get(), Spin::None, cd.get());
        ExpectMachineEqual(Vr, Vc, "DeltaFittedVxc");
    }
}


// The ASSEMBLY gate (Step 3c-2): the Hamiltonian's Ham_RealBlock fold over the term faces must equal the
// native complex assembly's real part BITWISE -- each term's block already does (the 3c-1 gates), and the
// fold accumulates them elementwise in the same order.  This is the exact matrix a tIrrepWF<double> child
// receives from its CalculateH inside a complex run.
namespace { struct MixedHam : Hamiltonian::tHamiltonianImp<dcmplx> {}; }

TEST(RealComplexTerms, HamiltonianAssemblyServesTheRealBlockBitwise)
{
    Rig rig;
    auto cd = rig.MakeDensity();
    const auto gth = Pseudopotential::GetGTH("Si","LDA",4);

    MixedHam ham;
    ham.Add(new Kinetic<dcmplx>);
    ham.Add(new Ven_PP_Short   (rig.st, &gth.local));
    ham.Add(new Ven_PP_Long    (rig.st, &gth.local));
    ham.Add(new Ven_PP_NonLocal(rig.st, &gth.nonlocal));
    ham.Add(new Vee_Hartree(Vee_Hartree::fbs_t(rig.cx->CreateCDFitBasisSet(rig.st.get(), qcMesh::MeshParams{}))));

    const chmat_t Hc = ham.GetMatrix(rig.cx.get(), Spin::None, cd.get(), nullptr);
    auto& rb = dynamic_cast<Hamiltonian::Ham_RealBlock&>(ham);   // the face a real WF child drives
    const rsmat_t Hr = rb.GetMatrix(rig.re.get(), Spin::None, cd.get(), nullptr);
    ExpectBitwiseEqual(Hr, Hc, "Ham_RealBlock fold");
}

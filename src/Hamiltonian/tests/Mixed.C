// File: src/Hamiltonian/tests/Mixed.C  THE LINEAR MIXER TERM (V-HAM-COMPOSE step 3, lane S 2026-11-26): sum_i w_i T_i over any terms.
//
// The gates the design note states (§3(c)): 1*T + 0*T' reproduces T BIT-IDENTICALLY (E and every block), and a*T + (1-a)*T = T at fixed D
// to roundoff.  On REAL terms: the PBE GGA term on water DZVP (GGA_Term.C's rig), two independent instances so nothing is shared.  Plus the
// plumbing that a real term cannot show cheaply: the family (an HF child makes the mixer whole-system), the one-mesh rule, the pin-28
// certificate surviving the scaling, and a child without a force being NAMED.
#include <cmath>
#include <cstdio>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <gtest/gtest.h>
import qchem.Hamiltonian.Internal.Mixed;             // the term under test (tests may cheat-import internals)
import qchem.Hamiltonian.Internal.Vxc_DFT_GGA;
import qchem.Hamiltonian.Internal.GGAFunctional;
import qchem.Hamiltonian.Internal.Libxc_GGA;
import qchem.Hamiltonian.Internal.Term;              // DensityMatrixRequirement
import qchem.Hamiltonian.Internal.GeometryGradient;  // tGeometryGradient_HT
import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.DensitySampler;
import qchem.ChargeDensity.Factory;
import qchem.BasisSet;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.Tests.PW_Fields;                        // OneBlockBS: a one-block whole basis for the refresh
import qchem.Mesh;
import qchem.Structure;
import qchem.Symmetry.Irrep;
import qchem.Types;
import qchem.Blaze;
import qchem.Structure.MeshFactory;
using namespace qchem;
using namespace qchem::Hamiltonian;

namespace
{
Molecule H2O()
{
    Molecule m;
    m.Insert(Site{8, rvec3_t( 0.000, 0.000, 0.000)});
    m.Insert(Site{1, rvec3_t( 1.430, 1.107, 0.000)});
    m.Insert(Site{1, rvec3_t(-1.430, 1.107, 0.000)});
    return m;
}
hmat_t<double> SomePSD_D(size_t n, double phase)
{
    hmat_t<double> D(n);
    for (size_t i=0; i<n; i++) for (size_t j=0; j<=i; j++)
    {
        double s=0.0;
        for (size_t m=0; m<5; m++) s+=std::sin(phase+1.1*i+0.7*m)*std::sin(phase+1.1*j+0.7*m);
        D(i,j)=2.0*s/double(n);
    }
    return D;
}

//! A whole-system term with a FIXED matrix M and energy E = 1/2 Tr(D M) under the name "Eee": the HF family's plumbing, no ERIs.
class StubHF : public virtual rDynamic_HF_HT
{
public:
    explicit StubHF(hmat_t<double> M) : itsM(std::move(M)) {}
    virtual const hmat_t<double>& GetMatrix(const robs_t*, const Spin&, DensityID, const rbs_t* whole) const override
    {
        if (!whole) throw std::logic_error("StubHF: no whole basis");
        return itsM;
    }
    virtual void RefreshForDensity(const rbs_t*, const rChargeDensity*) const override {++refreshes;}
    virtual void GetEnergy(EnergyBreakdown& e, const rDM_CD*) const override {e.Add("Eee", 1.25, EnergyRole::Potential, 2.5);}
    virtual std::ostream& Write(std::ostream& os) const override {return os << "StubHF";}
    mutable int refreshes=0;
private:
    hmat_t<double> itsM;
};

struct Rig
{
    Molecule                                         mol=H2O();
    std::unique_ptr<BasisSet::Real_BS>               bs;
    const BasisSet::Orbital_1E_IBS<double>*          block=nullptr;
    std::shared_ptr<const BasisSet::rFIT_SF_ABS>     fb;
    std::shared_ptr<const ChargeDensity::tDensityGradientSampler<double>> sampler;
    std::shared_ptr<const CompositeGGAFunctional>    pbe;
    Rig()
    {
        bs.reset(BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::DZVP, &mol,
                                             BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
        block=dynamic_cast<const BasisSet::Orbital_1E_IBS<double>*>((*bs)[0]);
        const qcMesh::MeshParams mp{.radial=qcMesh::RadialKind::MHL, .nRadial=40, .mhl_m=3, .mhl_alpha=2.0,
                                    .angular=qcMesh::AngularKind::Lebedev, .angularDegree=17, .beckeOrder=2};
        BasisSet::FitQuadrature quad=MakeMeshMaker(mol)->Make(mp, qcMesh::NoSymmetry());
        fb.reset(bs->CreateVxcFitBasisSet(&mol, BasisSet::VxcFit::Delta, &quad));
        sampler=ChargeDensity::MakeDensityGradientSampler(fb, quad);
        pbe=std::make_shared<const CompositeGGAFunctional>(std::vector<std::shared_ptr<const GGAFunctional>>{
            std::make_shared<const Libxc_GGA>(101), std::make_shared<const Libxc_GGA>(130)});   // GGA_X_PBE + GGA_C_PBE
    }
    rTerm GGA() const {return rTerm::Of(new tVxc_DFT_GGA<double>(pbe, sampler, SpinGroup::UnPolarized));}
    std::unique_ptr<ChargeDensity::tDM_CD<double>> Density(const hmat_t<double>& D) const
    {
        return std::unique_ptr<ChargeDensity::tDM_CD<double>>(ChargeDensity::IrrepCD_Factory<double>(D, block, bs->GetIrreps(Spin::None)[0], qchem::ChargeDensity::RhoRoute::PivotedCholesky));
    }
    //! A dynamic term's energy and its potential + energy blocks at \a cd, through the faces the Hamiltonian reads.
    struct Out { EnergyBreakdown e; hmat_t<double> V, E; };
    Out Evaluate(const tDynamic_HT<double>& t, const ChargeDensity::tDM_CD<double>* cd) const
    {
        qchem::tests::pw::OneBlockBS<double> one(block);
        t.RefreshForDensity(cd);
        t.BuildBlocks(&one, cd);
        Out o;
        t.GetEnergy(o.e, cd);
        o.V=t.GetMatrix (block, Spin::None, DensityID(cd));
        o.E=t.GetEMatrix(block, Spin::None, DensityID(cd));
        return o;
    }
};
double MaxAbsDiff(const hmat_t<double>& a, const hmat_t<double>& b)
{
    double d=0.0;
    for (size_t i=0; i<a.rows(); i++) for (size_t j=0; j<a.columns(); j++) d=std::max(d, std::fabs(a(i,j)-b(i,j)));
    return d;
}
} //anon

// 1*T + 0*T' = T, BIT FOR BIT: 1.0*x is x, and adding 0.0*y (y finite) is adding +0 -- so nothing in the mixer may reorder or round.
TEST(Mixed, OneTimesATermPlusZeroTimesAnotherIsTheTermBitForBit)
{
    Rig rig;
    const auto cd=rig.Density(SomePSD_D(rig.block->GetNumFunctions(), 0.3));
    rTerm plain=rig.GGA();
    std::vector<tWeightedTerm<double>> parts;
    parts.push_back({1.0, rig.GGA()});
    parts.push_back({0.0, rig.GGA()});
    rTerm mixed=Mix(std::move(parts));
    ASSERT_TRUE(mixed.dynamic) << "two dynamic children make a dynamic mixer";
    const auto a=rig.Evaluate(*plain.dynamic, cd.get()), b=rig.Evaluate(*mixed.dynamic, cd.get());
    EXPECT_EQ(a.e["Exc"], b.e["Exc"]);
    EXPECT_EQ(MaxAbsDiff(a.V, b.V), 0.0) << "the potential block";
    EXPECT_EQ(MaxAbsDiff(a.E, b.E), 0.0) << "the energy block";
    EXPECT_EQ(a.e.charge.xcRhoCertified, b.e.charge.xcRhoCertified) << "the pin-28 certificate survives the scaling";
    EXPECT_TRUE(b.e.charge.xcBuilt);
}

// a*T + (1-a)*T = T at fixed D, to roundoff (the weights sum to one; the two children are independent instances).
TEST(Mixed, AWeightedSumOfATermWithItselfIsTheTerm)
{
    Rig rig;
    const auto cd=rig.Density(SomePSD_D(rig.block->GetNumFunctions(), 0.9));
    rTerm plain=rig.GGA();
    const double a=0.3;
    std::vector<tWeightedTerm<double>> parts;
    parts.push_back({a,     rig.GGA()});
    parts.push_back({1.0-a, rig.GGA()});
    rTerm mixed=Mix(std::move(parts));
    const auto x=rig.Evaluate(*plain.dynamic, cd.get()), y=rig.Evaluate(*mixed.dynamic, cd.get());
    std::printf("  [mixer] E_xc %.14f vs %.14f, max|dV| %.2e, max|dE| %.2e\n", x.e["Exc"], y.e["Exc"], MaxAbsDiff(x.V,y.V), MaxAbsDiff(x.E,y.E));
    EXPECT_NEAR(y.e["Exc"], x.e["Exc"], 1e-14*std::fabs(x.e["Exc"]));
    EXPECT_LT(MaxAbsDiff(x.V, y.V), 1e-14);
    EXPECT_LT(MaxAbsDiff(x.E, y.E), 1e-14);
    std::ostringstream os; os << *mixed.dynamic;
    EXPECT_NE(os.str().find("Linear mix of 2 terms"), std::string::npos) << os.str();
}

// THE FAMILY: an HF child makes the mixer whole-system; its block is the weighted sum (the dynamic child's block + w * the HF matrix);
// its energy entries are the children's, scaled; each child is refreshed once per refresh.
TEST(Mixed, AnHFChildMakesTheMixerWholeSystemAndItsBlockIsTheWeightedSum)
{
    Rig rig;
    const size_t n=rig.block->GetNumFunctions();
    const auto cd=rig.Density(SomePSD_D(n, 0.5));
    hmat_t<double> M=SomePSD_D(n, 2.1);
    auto* stub=new StubHF(M);
    rTerm plain=rig.GGA();
    std::vector<tWeightedTerm<double>> parts;
    parts.push_back({1.0, rig.GGA()});
    parts.push_back({0.25, rTerm::Of(stub)});
    rTerm mixed=Mix(std::move(parts));
    ASSERT_TRUE(mixed.hf) << "an HF child makes the mixer whole-system";
    ASSERT_FALSE(mixed.dynamic);
    const auto ref=rig.Evaluate(*plain.dynamic, cd.get());
    qchem::tests::pw::OneBlockBS<double> one(rig.block);
    mixed.hf->RefreshForDensity(&one, cd.get());
    EXPECT_EQ(stub->refreshes, 1);
    const hmat_t<double>& H=mixed.hf->GetMatrix(rig.block, Spin::None, DensityID(cd.get()), &one);
    hmat_t<double> want=ref.V; want+=0.25*M;
    EXPECT_LT(MaxAbsDiff(H, want), 1e-15);
    EnergyBreakdown e;
    mixed.hf->GetEnergy(e, cd.get());
    EXPECT_EQ(e["Eee"], 0.25*1.25);
    EXPECT_EQ(e["Exc"], ref.e["Exc"]);
    const auto* r=dynamic_cast<const DensityMatrixRequirement*>(mixed.hf.get());
    ASSERT_TRUE(r);
    EXPECT_TRUE(r->RequiresDensityMatrix()) << "an HF child that does not answer needs D (tHamiltonianImp::Add(HF)'s default)";
    EXPECT_THROW((void)mixed.hf->GetMatrix(rig.block, Spin::None, DensityID(rig.Density(SomePSD_D(n,0.7)).get()), &one), std::logic_error)
        << "a block read for a density the mixer was not refreshed for";
}

// A child with no force face is NAMED by MissingGradient (never a silent zero); refusals at construction.
TEST(Mixed, AChildWithoutAForceIsNamedAndBadMixesAreRefused)
{
    Rig rig;
    std::vector<tWeightedTerm<double>> parts;
    parts.push_back({1.0, rig.GGA()});
    parts.push_back({0.5, rTerm::Of(new StubHF(SomePSD_D(rig.block->GetNumFunctions(), 1.0)))});
    rTerm mixed=Mix(std::move(parts));
    const auto* g=dynamic_cast<const tGeometryGradient_HT<double>*>(mixed.hf.get());
    ASSERT_TRUE(g);
    const std::string why=g->MissingGradient(rig.bs.get(), &rig.mol);
    EXPECT_FALSE(why.empty());
    std::printf("  [mixer] MissingGradient: %s\n", why.c_str());

    EXPECT_THROW(Mix(std::vector<tWeightedTerm<double>>{}), std::invalid_argument);
    std::vector<tWeightedTerm<double>> nan;
    nan.push_back({std::nan(""), rig.GGA()});
    EXPECT_THROW(Mix(std::move(nan)), std::invalid_argument);
    std::vector<tWeightedTerm<double>> twoMeshes;
    rTerm a=rig.GGA(), b=rig.GGA();
    a.mesh=rig.sampler; b.mesh=rig.sampler;
    twoMeshes.push_back({1.0, std::move(a)});
    twoMeshes.push_back({1.0, std::move(b)});
    EXPECT_THROW(Mix(std::move(twoMeshes)), std::logic_error) << "one Hamiltonian, one mesh";
}

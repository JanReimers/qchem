// File: Response/tests/Kernel.C  The response kernel against its FINITE-DIFFERENCE oracle
// (doc/LinearResponsePlan.md §3 H1/H2, ruling D6).
//
// THE ORACLE.  δF = [F(D0 + hδD) - F(D0 - hδD)] / 2h, built from the Hamiltonian's PUBLIC GetMatrix and
// ordinary densities: it knows nothing about how any term linearises, so it checks every analytic
// tResponse_HT the same way.  It lives here, in the test tree (D6), and reads δD's AO blocks through the one
// friend src/forward.H grants it -- the TransitionDensity face never hands its matrices out.
//
// For HF the gate needs NO SCF: F = h + J[D] - K[D] is exactly linear in D, so the central difference is exact
// (to rounding) at ANY D0.  Random D0 and random δD therefore test the whole chain -- transition density ->
// HF sweep face -> the term's response slot -> the kernel fold -> the TransitionFock -- on real ERIs.
#include "gtest/gtest.h"
#include <cmath>
#include <memory>
#include <random>
#include <stdexcept>
#include <vector>
#include <nlohmann/json.hpp>
#include "forward.H"
import qchem.Hamiltonian.Factory;
import qchem.BasisSet.Gaussian.Point.Factory;
import qchem.ChargeDensity.Factory;                   // IrrepCD_Factory, RhoRoute
import qchem.CompositeCD;
import qchem.ChargeDensity.Internal.TransitionDensity; // the concrete the friend reads (tests may cheat)
import qchem.Structure;
import qchem.Mesh;                                    // qcMesh::MeshParams (the LDA Hamiltonian's XC mesh)
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using ChargeDensity::TransitionBlock;
using ChargeDensity::TransitionDensity;

//! THE FRIEND (src/forward.H): the one door to a transition density's AO blocks.
class TransitionDensityTests
{
public:
    template <class T> static const std::vector<TransitionBlock<T>>& Blocks(const TransitionDensity<T>& d)
    {
        auto* c=dynamic_cast<const ChargeDensity::AO_TransitionDensity<T>*>(&d);   // tests may cheat: abstract -> concrete
        if (!c) throw std::logic_error("TransitionDensityTests: not an AO transition density");
        if (c->itsBlocks.empty()) throw std::logic_error("TransitionDensityTests: a channel VIEW has no blocks of its own");
        return c->itsBlocks;
    }
};

namespace {

//! The finite-difference response kernel (D6): the oracle for every analytic tResponse_HT.
template <class T> class FD_ResponseKernel : public Hamiltonian::ResponseKernel<T>
{
public:
    //! \a D0 = the linearisation point, as AO blocks in the SAME order a transition density will list them.
    FD_ResponseKernel(Hamiltonian::tHamiltonian<T>& H, const Hamiltonian::tbs_t<T>* wholeBasis,
                      std::vector<TransitionBlock<T>> D0, double h)
        : itsH(&H), itsWB(wholeBasis), itsD0(std::move(D0)), itsH_(h) {}

    virtual std::unique_ptr<Hamiltonian::TransitionFock<T>> InducedFock(const TransitionDensity<T>& delta) const override
    {
        const auto& dD=TransitionDensityTests::Blocks(delta);
        if (dD.size()!=itsD0.size()) throw std::logic_error("FD_ResponseKernel: δD and D0 have different block lists");
        auto Dp=Displaced(dD, +itsH_), Dm=Displaced(dD, -itsH_);
        auto dF=std::make_unique<Hamiltonian::AO_TransitionFock<T>>(delta.Coupling());
        for (const auto& b : itsD0)
        {
            const hmat_t<T> Fp=itsH->GetMatrix(b.bs, b.irrep.ms, Dp.get(), itsWB);
            const hmat_t<T> Fm=itsH->GetMatrix(b.bs, b.irrep.ms, Dm.get(), itsWB);
            dF->Add(b.irrep, b.irrep, mat_t<T>((Fp-Fm)/T(2.0*itsH_)));
        }
        return dF;
    }
private:
    std::unique_ptr<ChargeDensity::tComposite_CD<T>> Displaced(const std::vector<TransitionBlock<T>>& dD, double h) const
    {
        auto cd=std::make_unique<ChargeDensity::tComposite_CD<T>>();
        for (size_t b=0;b<itsD0.size();b++)
        {
            if (dD[b].irrep<itsD0[b].irrep || itsD0[b].irrep<dD[b].irrep)
                throw std::logic_error("FD_ResponseKernel: δD and D0 list their blocks in a different order");
            const hmat_t<T> D=itsD0[b].dD+T(h)*dD[b].dD;
            cd->Insert(std::unique_ptr<ChargeDensity::tDM_CD<T>>(
                ChargeDensity::IrrepCD_Factory<T>(D, itsD0[b].bs, itsD0[b].irrep, ChargeDensity::RhoRoute::Direct)), itsD0[b].irrep);
        }
        return cd;
    }
    Hamiltonian::tHamiltonian<T>*       itsH;   // GetMatrix is non-const on the face
    const Hamiltonian::tbs_t<T>*        itsWB;
    std::vector<TransitionBlock<T>>     itsD0;
    double                              itsH_;
};

std::shared_ptr<const Structure> Water()     // M_Calculation's geometry, bohr
{
    auto w=std::make_shared<Molecule>();
    w->Insert(new Atom(8, 0, Vector3D<double>(0,  0.0,   0.0)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0,  1.431, 1.107)));
    w->Insert(new Atom(1, 0, Vector3D<double>(0, -1.431, 1.107)));
    return w;
}

rsmat_t RandomSymmetric(size_t n, std::mt19937& g, double scale)
{
    std::uniform_real_distribution<double> u(-scale, scale);
    rsmat_t M(n);
    for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) M(i,j)=u(g);
    return M;
}

//! Build one block list over the whole basis for the spin irreps of \a g, filled by \a fill.
template <class Fill> std::vector<TransitionBlock<double>> Blocks(const BasisSet::Real_BS& bs, SpinGroup g, Fill&& fill)
{
    std::vector<TransitionBlock<double>> out;
    const std::vector<Spin> spins = g==SpinGroup::Polarized ? std::vector<Spin>{Spin::Up, Spin::Down}
                                                            : std::vector<Spin>{Spin::None};
    for (const Spin& s : spins)
        for (const auto* ob : bs.Iterate<Hamiltonian::robs_t>())
            out.push_back({ob->GetIrrep(s), ob, fill(ob->GetNumFunctions())});
    return out;
}

void AnalyticHF_EqualsFD(SpinGroup g)
{
    auto st=Water();
    std::unique_ptr<BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(
        nlohmann::json{{"basis","dzvp"},{"engine","mnd"},{"angular","cartesian"}}, st.get()));
    std::unique_ptr<Hamiltonian::rHamiltonian> H(Hamiltonian::Factory(Hamiltonian::Model::HF, g, st));
    std::mt19937 rng(20260928);
    auto D0=Blocks(*bs, g, [&](size_t n){return RandomSymmetric(n, rng, 0.3);});
    auto dD=Blocks(*bs, g, [&](size_t n){return RandomSymmetric(n, rng, 1.0);});

    auto D0cd=std::make_unique<ChargeDensity::rComposite_CD>();
    for (const auto& b : D0)
        D0cd->Insert(std::unique_ptr<ChargeDensity::rDM_CD>(
            ChargeDensity::IrrepCD_Factory<double>(b.dD, b.bs, b.irrep, ChargeDensity::RhoRoute::Direct)), b.irrep);

    auto rule=std::make_shared<Symmetry::Invariant>();
    auto delta=ChargeDensity::AO_TransitionDensity_Factory<double>(dD, rule);
    auto analytic=H->MakeResponseKernel(bs.get(), D0cd.get());
    FD_ResponseKernel<double> fd(*H, bs.get(), D0, 1e-3);
    auto A=analytic->InducedFock(*delta);
    auto F=fd.InducedFock(*delta);

    for (const auto& b : dD)
    {
        const rmat_t a=A->Matrix(b.irrep, b.irrep), f=F->Matrix(b.irrep, b.irrep);
        double diff=0, scale=0;
        for (size_t i=0;i<f.rows();i++)
            for (size_t j=0;j<f.columns();j++) {diff=std::max(diff, std::fabs(a(i,j)-f(i,j))); scale=std::max(scale, std::fabs(f(i,j)));}
        EXPECT_GT(scale, 1e-2) << b.irrep;      // not a vacuous comparison of two zeros
        EXPECT_LT(diff/scale, 1e-9) << b.irrep << ": analytic J/K kernel vs finite difference";
    }
}

} // namespace

TEST(ResponseKernel, HF_AnalyticEqualsFiniteDifference_UnPol) {AnalyticHF_EqualsFD(SpinGroup::UnPolarized);}
TEST(ResponseKernel, HF_AnalyticEqualsFiniteDifference_Pol)   {AnalyticHF_EqualsFD(SpinGroup::Polarized);}

//! A dynamic term with no response capability makes MakeResponseKernel THROW, naming it (ruling Q2) -- the
//! molecular LDA Hamiltonian until R2 gives its fitted Coulomb/XC terms the face.
TEST(ResponseKernel, MissingCapabilityThrowsAtConstruction)
{
    auto st=Water();
    std::unique_ptr<BasisSet::Real_BS> bs(BasisSet::Gaussian::Factory(
        nlohmann::json{{"basis","dzvp"},{"engine","mnd"},{"angular","cartesian"}}, st.get()));
    std::unique_ptr<Hamiltonian::rHamiltonian> H(Hamiltonian::Factory(Hamiltonian::Model::LDA, SpinGroup::UnPolarized, st,
        qcMesh::MeshParams{.radial=qcMesh::RadialKind::MHL, .nRadial=30, .mhl_m=3, .mhl_alpha=2.0,
                           .angular=qcMesh::AngularKind::Lebedev, .angularDegree=5, .beckeOrder=2}, bs.get(), 0.7));
    try
    {
        (void)H->MakeResponseKernel(bs.get(), nullptr);
        FAIL() << "an LDA Hamiltonian has no response kernel yet: MakeResponseKernel must throw";
    }
    catch (const std::logic_error& e)
    {
        EXPECT_NE(std::string(e.what()).find("cannot be linearised"), std::string::npos) << e.what();
    }
}

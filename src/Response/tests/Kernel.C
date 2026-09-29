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
#include <cstdlib>
#include <complex>
#include <iostream>
#include <limits>
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
import qchem.Calculation;
import qchem.SolidCalculation;                        // the converged GPW solid the periodic gates run on (R2)
import qchem.Lattice_3D;
import qchem.BasisSet.Lattice.BasisSet;                             // the converged H2O the FD-driven polarisability runs on
import qchem.Response;                                // Reference / frame / probe / LinearResponse
import qchem.Orbitals;                                // the ground-state density blocks the FD kernel displaces
import qchem.SCFParams;
import qchem.Mesh;                                    // qcMesh::MeshParams (the LDA Hamiltonian's XC mesh)
import qchem.Types;
import qchem.Blaze;

using namespace qchem;
using ChargeDensity::TransitionBlock;
using ChargeDensity::TransitionDensity;

//! THE FRIEND of qchem::Calculation (src/forward.H): the converged run's Hamiltonian, wave function and basis.
class ResponseFacadeTests
{
public:
    static Hamiltonian::rHamiltonian&                   Ham  (Calculation& c)       {return *c.itsHam;}
    static const WaveFunction::tWaveFunction<double>&   WF   (const Calculation& c) {return *c.itsScf->GetWaveFunction();}
    static const BasisSet::Real_BS*                     Basis(const Calculation& c) {return c.itsBasis;}
    static const Structure&                             St   (const Calculation& c) {return *c.itsStructure;}
    // The periodic facade (R2): its Hamiltonian and wave function, through its private door.
    static Hamiltonian::cHamiltonian&                   Ham  (const SolidCalculation& c) {return c.ResponseHamiltonian();}
    static const WaveFunction::cWaveFunction&           WF   (const SolidCalculation& c) {return c.ResponseWaveFunction();}
};

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

//! One block of a GROUND-STATE-shaped density (Hermitian D on one block): the FD oracle's linearisation point and
//! the q = 0 δD's these gates draw.  \c AsTransition makes the q = 0 transition density of such a list.
template <class T> struct GroundBlock
{
    Irrep                           irrep;
    const ChargeDensity::tobs_t<T>* bs=nullptr;
    hmat_t<T>                       D;
};
template <class T> std::vector<TransitionBlock<T>> AsTransition(const std::vector<GroundBlock<T>>& g)
{
    std::vector<TransitionBlock<T>> out;
    for (const auto& b : g) out.push_back({b.irrep, b.irrep, b.bs, b.bs, mat_t<T>(b.D)});
    return out;
}

//! The finite-difference response kernel (D6): the oracle for every analytic tResponse_HT.
template <class T> class FD_ResponseKernel : public Hamiltonian::ResponseKernel<T>
{
public:
    //! \a D0 = the linearisation point, as AO blocks in the SAME order a transition density will list them.
    FD_ResponseKernel(Hamiltonian::tHamiltonian<T>& H, const Hamiltonian::tbs_t<T>* wholeBasis,
                      std::vector<GroundBlock<T>> D0, double h)
        : itsH(&H), itsWB(wholeBasis), itsD0(std::move(D0)), itsH_(h) {}

    virtual std::unique_ptr<Hamiltonian::TransitionFock<T>> InducedFock(const TransitionDensity<T>& delta) const override
    {
        const auto& dD=TransitionDensityTests::Blocks(delta);
        if (dD.size()!=itsD0.size()) throw std::logic_error("FD_ResponseKernel: δD and D0 have different block lists");
        // The FOUR-POINT stencil, O(h^4): [8(F(h)-F(-h)) - (F(2h)-F(-2h))] / 12h.  Measured on GPW Si (R2): the
        // two-point form bottoms out at ~2e-6 relative (h^2 truncation meets rounding near h = 1e-4), which is an
        // ORACLE floor, not a kernel error; four points move the floor below 1e-7.  (Exact for HF either way.)
        auto dF=std::make_unique<Hamiltonian::AO_TransitionFock<T>>(delta.Coupling());
        auto D1p=Displaced(dD, +itsH_), D1m=Displaced(dD, -itsH_), D2p=Displaced(dD, +2*itsH_), D2m=Displaced(dD, -2*itsH_);
        for (const auto& b : itsD0)
        {
            auto F=[&](const auto& D){return mat_t<T>(itsH->GetMatrix(b.bs, b.irrep.ms, D.get(), itsWB));};
            const mat_t<T> d1=F(D1p)-F(D1m), d2=F(D2p)-F(D2m);
            dF->Add(b.irrep, b.irrep, mat_t<T>((T(8.0)*d1-d2)/T(12.0*itsH_)));
        }
        return dF;
    }
private:
    std::unique_ptr<ChargeDensity::tComposite_CD<T>> Displaced(const std::vector<TransitionBlock<T>>& dD, double h) const
    {
        auto cd=std::make_unique<ChargeDensity::tComposite_CD<T>>();
        for (size_t b=0;b<itsD0.size();b++)
        {
            if (dD[b].ket<itsD0[b].irrep || itsD0[b].irrep<dD[b].ket)
                throw std::logic_error("FD_ResponseKernel: δD and D0 list their blocks in a different order");
            if (dD[b].bra.SequenceIndex()!=dD[b].ket.SequenceIndex())
                throw std::logic_error("FD_ResponseKernel: a (bra != ket) pair -- the FD oracle is q = 0 only");
            hmat_t<T> D=itsD0[b].D;
            for (size_t i=0;i<D.rows();i++) for (size_t j=i;j<D.columns();j++) D(i,j)+=T(h)*dD[b].dD(i,j);
            cd->Insert(std::unique_ptr<ChargeDensity::tDM_CD<T>>(
                ChargeDensity::IrrepCD_Factory<T>(D, itsD0[b].bs, itsD0[b].irrep, ChargeDensity::RhoRoute::Direct)), itsD0[b].irrep);
        }
        return cd;
    }
    Hamiltonian::tHamiltonian<T>*       itsH;   // GetMatrix is non-const on the face
    const Hamiltonian::tbs_t<T>*        itsWB;
    std::vector<GroundBlock<T>>         itsD0;
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
template <class Fill> std::vector<GroundBlock<double>> Blocks(const BasisSet::Real_BS& bs, SpinGroup g, Fill&& fill)
{
    std::vector<GroundBlock<double>> out;
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
            ChargeDensity::IrrepCD_Factory<double>(b.D, b.bs, b.irrep, ChargeDensity::RhoRoute::Direct)), b.irrep);

    auto rule=std::make_shared<Symmetry::Invariant>();
    auto delta=ChargeDensity::AO_TransitionDensity_Factory<double>(AsTransition(dD), rule);
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

//=====================================================================================================
//  THE SOLVER WITH THE ORACLE KERNEL.  The same Reference / frame / dipole probe / GMRES as
//  Calculation::StaticPolarizability, but the kernel is the finite difference of the Hamiltonian's own
//  GetMatrix about the CONVERGED density.  Two uses:
//   * HF: it must reproduce the analytic kernel's polarisability -- the FD kernel is trusted INSIDE the solver;
//   * LDA: the one route to an LDA polarisability until R2 gives the fitted terms an analytic face.
//=====================================================================================================
namespace {

const SCFParams tight = {.NMaxIter=80, .MinΔρ=1e-9, .MinΔFD=1e-10, .MinVirial=1e2};
const qcMesh::MeshParams dipoleMesh={.radial=qcMesh::RadialKind::MHL, .nRadial=80, .mhl_m=3, .mhl_alpha=2.0,
                                     .angular=qcMesh::AngularKind::Lebedev, .angularDegree=35, .beckeOrder=3};

//! The converged density matrix of every block, \f$\sum_i n_i c_ic_i^\dagger\f$, in the wave function's order.
std::vector<GroundBlock<double>> GroundDensity(const WaveFunction::tWaveFunction<double>& wf)
{
    std::vector<GroundBlock<double>> out;
    for (const Irrep& ir : wf.GetQNs())
    {
        const auto* os=dynamic_cast<const Orbitals::TOrbitals<double>*>(wf.GetOrbitals(ir));
        const auto* bs=dynamic_cast<const Hamiltonian::robs_t*>(os->GetBasisSet());
        const size_t n=bs->GetNumFunctions();
        rsmat_t D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=0.0;
        for (const auto* o : os->Iterate<Orbitals::TOrbital<double>>())
        {
            const double occ=o->GetOccupation();
            if (occ==0.0) continue;
            const rvec_t& c=o->GetCoeff();
            for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)+=occ*c[i]*c[j];
        }
        out.push_back({ir, bs, D});
    }
    return out;
}

//! alpha = -chi of the dipole channels, with the finite-difference kernel of step \a h.
rmat_t FD_Polarizability(Calculation& calc, double h)
{
    const auto& wf=ResponseFacadeTests::WF(calc);
    const Response::Reference ref=Response::MakeReference(wf, OccupationConfig{}, {.acrossK=true, .acrossSpin=false},
                                                          std::numeric_limits<double>::quiet_NaN());
    const auto frame=Response::MakeOrbitalFrame(ref, wf);
    const auto probe=Response::MakeDipoleProbe(ref, frame, wf, ResponseFacadeTests::St(calc).CreateIntegrationMesh(dipoleMesh));
    FD_ResponseKernel<double> fd(ResponseFacadeTests::Ham(calc), ResponseFacadeTests::Basis(calc), GroundDensity(wf), h);
    auto r=Response::LinearResponse(ref, frame, fd, probe, std::make_shared<Symmetry::Invariant>(), {.tol=1e-9});
    EXPECT_TRUE(r.IsOk()) << (r ? "" : r.Error().detail);
    rmat_t a(3,3,0.0);
    if (!r) return a;
    r->Write(std::cout);
    for (size_t i=0;i<3;i++) for (size_t j=0;j<3;j++) a(i,j)=-r->chi(i,j).real();
    return a;
}

} // namespace

//! The FD kernel INSIDE the solver reproduces the analytic HF polarisability (and so PySCF's).
TEST(ResponsePolarizability, HF_FD_Kernel_eqAnalytic)
{
    Calculation calc(*Water(), {.basis="dzvp"});
    ASSERT_TRUE(calc.Converge(tight));
    auto analytic=calc.StaticPolarizability(dipoleMesh);
    ASSERT_TRUE(analytic.IsOk()) << analytic.Error().detail;
    const rmat_t fd=FD_Polarizability(calc, 1e-3);
    for (size_t i=0;i<3;i++) EXPECT_NEAR(fd(i,i), analytic.Value()(i,i), 1e-7*analytic.Value()(i,i)) << "alpha_" << i << i;
}

//! LDA (Slater + VWN5) CPKS through the FD kernel -- the only LDA kernel until R2.  A LOOSE oracle, and the
//! looseness is MEASURED to be the GROUND STATE's XC route, not the response (2026-09-28, H2O/dzvp):
//!     XC mesh                    E (Ha)       alpha xx / yy / zz  vs PySCF lda,vwn5 CPKS
//!     default (MHL 30 / Leb 5)   -75.93246    -0.2%  -1.2%  -2.0%
//!     MHL 80 / Lebedev 35        -75.8726     +0.6%  +1.4%  +0.7%      (31 s -- too slow for the gate)
//!     PySCF                      -75.87730    (3.44581  7.33076  5.91599)
//! The residual ~1% at a fine mesh is our FITTED Coulomb/XC against PySCF's exact J.  The FD step is NOT a
//! factor: h = 1e-3 and 1e-4 agree to 1e-6 (checked below).  Gated at 3% on the default mesh.
TEST(ResponsePolarizability, LDA_FD_Kernel_vsPySCF)
{
    const double pyscf[3]={3.44581010, 7.33075620, 5.91599197};   // PySCF RKS lda,vwn / dzvp(cart) CPKS
    Calculation calc(*Water(), {.basis="dzvp", .model=Hamiltonian::Model::LDA});
    ASSERT_TRUE(calc.Converge(tight));
    const rmat_t a=FD_Polarizability(calc, 1e-4);
    const rmat_t b=FD_Polarizability(calc, 1e-3);
    for (size_t i=0;i<3;i++)
    {
        EXPECT_NEAR(a(i,i), b(i,i), 1e-5*a(i,i)) << "alpha_" << i << i << " depends on the FD step: the kernel is not linearised";
        std::cout << "[LDA alpha] " << i << i << "  ours " << a(i,i) << "  PySCF " << pyscf[i]
                  << "  rel " << (a(i,i)-pyscf[i])/pyscf[i] << std::endl;
        EXPECT_NEAR(a(i,i), pyscf[i], 0.03*pyscf[i]) << "alpha_" << i << i;
    }
}

//=====================================================================================================
//  R2: THE PERIODIC KERNEL (GPW Hartree + ALDA f_xc, +U at U = 0) against the same FD oracle, on a CONVERGED
//  solid -- XC is nonlinear, so D0 must be a real, positive density (random D0 would sample f_xc at rho < 0).
//  Si diamond at Gamma (SCFTrace's cell), complex ansatz (the response faces serve complex blocks), and a
//  Si-p manifold at U = 0 so the run carries Hubbard channels for the facade's probe.
//=====================================================================================================
namespace {

//! \a U (Hartree) and \a alpha on site 0's p manifold; \a spectator adds site 1's p at U = 0 (a projector-set
//! member that carries no +U -- the §3d Q10 case: measured, not perturbed by default).
SCFParams SiParams()
{
    SCFParams par;
    // 200, not 80: the POLARIZED k211 ground state needs more than 80 (2026-09-29); a cap costs nothing unless hit (trap 3).
    par.NMaxIter=200; par.MinΔρ=1e-7; par.MinΔE=1e-10; par.MinΔFD=1e-7; par.MinVirial=1e30; par.MinFD=1e30;
    return par;
}
std::unique_ptr<SolidCalculation> ConvergedSi(SpinGroup g, double alpha=0.0, double U=0.0, bool spectator=false,
                                              ivec3_t kmesh=ivec3_t(1,1,1), bool impose=false)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, kmesh);
    auto mol=std::shared_ptr<const BasisSet::Real_BS>(
        BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    const SCFParams par=SiParams();
    SolidCalcOptions o{.Nelec=8, .multiplicity = g==SpinGroup::Polarized ? 1 : 0, .species={{"Si",4}},
                       .densityEcut=20.0, .hubbard={{.site=0, .l=1, .U=U, .alpha=alpha}}, .forceComplex=true};
    if (spectator) o.hubbard.push_back({.site=1, .l=1, .U=0.0});
    o.imposeSymmetry=impose;
    auto calc=std::make_unique<SolidCalculation>(lat, mol, o, par);
    EXPECT_TRUE(calc->DidConverge());
    return calc;
}

//! The converged density matrix of every (complex) block, in the wave function's order.
std::vector<GroundBlock<dcmplx>> GroundDensity(const WaveFunction::cWaveFunction& wf)
{
    std::vector<GroundBlock<dcmplx>> out;
    for (const Irrep& ir : wf.GetQNs())
    {
        const auto* os=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(wf.GetOrbitals(ir));
        if (!os) throw std::logic_error("GroundDensity: a real block -- run the ground state with forceComplex");
        const auto* bs=dynamic_cast<const Hamiltonian::cobs_t*>(os->GetBasisSet());
        const size_t n=bs->GetNumFunctions();
        chmat_t D(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)=0.0;
        for (const auto* o : os->Iterate<Orbitals::TOrbital<dcmplx>>())
        {
            const double occ=o->GetOccupation();
            if (occ==0.0) continue;
            const cvec_t& c=o->GetCoeff();
            for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) D(i,j)+=occ*c[i]*std::conj(c[j]);
        }
        out.push_back({ir, bs, D});
    }
    return out;
}

void PeriodicAnalyticEqualsFD(SpinGroup g)
{
    auto calc=ConvergedSi(g);
    const auto& wf=ResponseFacadeTests::WF(*calc);
    auto& H=ResponseFacadeTests::Ham(*calc);
    const auto D0=GroundDensity(wf);
    std::mt19937 rng(20260928);
    std::uniform_real_distribution<double> u(-0.05, 0.05);
    std::vector<GroundBlock<dcmplx>> dD;
    for (const auto& b : D0)
    {
        const size_t n=b.bs->GetNumFunctions();
        chmat_t X(n);
        for (size_t i=0;i<n;i++) {X(i,i)=u(rng); for (size_t j=i+1;j<n;j++) X(i,j)=dcmplx(u(rng),u(rng));}
        dD.push_back({b.irrep, b.bs, X});
    }
    const auto D0cd=wf.GetChargeDensity();
    auto analytic=H.MakeResponseKernel(&calc->Basis(), D0cd.get());
    FD_ResponseKernel<dcmplx> fd(H, &calc->Basis(), D0, 1e-3);   // 4-point stencil: h=1e-3 is its sweet spot (measured)
    auto delta=ChargeDensity::AO_TransitionDensity_Factory<dcmplx>(AsTransition(dD), std::make_shared<Symmetry::Invariant>());
    auto A=analytic->InducedFock(*delta);
    auto F=fd.InducedFock(*delta);
    for (const auto& b : dD)
    {
        const mat_t<dcmplx> a=A->Matrix(b.irrep, b.irrep), f=F->Matrix(b.irrep, b.irrep);
        double diff=0, scale=0;
        for (size_t i=0;i<f.rows();i++)
            for (size_t j=0;j<f.columns();j++) {diff=std::max(diff, std::abs(a(i,j)-f(i,j))); scale=std::max(scale, std::abs(f(i,j)));}
        std::cout << "[R2 kernel] " << b.irrep << "  max|FD| " << scale << "  max|analytic-FD| " << diff
                  << "  rel " << diff/scale << std::endl;
        EXPECT_GT(scale, 1e-3) << b.irrep;
        // MEASURED 2026-09-28 (4-point FD, h=1e-3, uniform XC raster): UnPol 7.2e-7; Pol ~4e-6 in the perturbed-
        // spin block -- a floor flat in h (1e-4..1e-3), in the dD amplitude, and in every knob tried (+U off, D0
        // route, D-aware screen, screen eps 1e-14); the functional's own kernel passes XCKernel.* to 1e-7.  OPEN,
        // doc/LinearResponsePlan.md §5d; 1000x below what chi or U can resolve.  ⚠ Do NOT move this gate to the
        // Becke mesh with a RANDOM dD: its far-tail points put h*drho >> rho0 and the FD saturates on the rho>0
        // guards (measured 5e-3 at h*amp=5e-5, 3e-5 at 5e-7) -- an ORACLE limit, not a kernel error.
        EXPECT_LT(diff/scale, 1e-5) << b.irrep << ": analytic Hartree + ALDA kernel vs finite difference";
    }
}

} // namespace

TEST(ResponseKernel, GPW_Si_AnalyticEqualsFiniteDifference_UnPol) {PeriodicAnalyticEqualsFD(SpinGroup::UnPolarized);}
TEST(ResponseKernel, GPW_Si_AnalyticEqualsFiniteDifference_Pol)   {PeriodicAnalyticEqualsFD(SpinGroup::Polarized);}

//=====================================================================================================
//  THE KERNEL IS LINEAR -- at every SCALE of δD (LinearResponsePlan §3e open item).  GMRES assumes a fixed linear
//  map, and on NiO its recurrence estimate ran 100x below the TRUE residual.  Suspect: a screen with an ABSOLUTE
//  tolerance (the D-aware collocation screen), which a Krylov solve probes with unit-norm vectors of tiny
//  components.  The R2 FD gate varied the knobs but never the SCALE.  So: K[s δ]/s must equal K[δ] down to
//  s = 1e-8, and K[δ1+δ2] = K[δ1]+K[δ2].  Multi-k (w_k = 1/2), NiO's shape.
//=====================================================================================================
TEST(ResponseKernel, GPW_Si_k211_KernelIsLinear_AtEveryScale)
{
    // The U = 2 eV state, +U FROZEN: on SiParams the k211 Si ground state converges at U = 2 eV (59 iterations) but
    // not at U = 0 nor polarized (200 iterations, 2026-09-29 -- an SCF-recipe matter).  The frozen +U kernel is zero,
    // and the screen under test is spin-agnostic, so the claim is unchanged.
    auto calc=ConvergedSi(SpinGroup::UnPolarized, 0.0, 2.0/27.211386245988, false, ivec3_t(2,1,1));
    const auto& wf=ResponseFacadeTests::WF(*calc);
    auto& H=ResponseFacadeTests::Ham(*calc);
    H.GetHubbardUTarget()->FreezeOccupations(true);
    const auto D0=GroundDensity(wf);
    const auto D0cd=wf.GetChargeDensity();
    auto K=H.MakeResponseKernel(&calc->Basis(), D0cd.get());
    std::mt19937 rng(20260929);
    std::uniform_real_distribution<double> u(-0.05, 0.05);
    auto random=[&]()
    {
        std::vector<GroundBlock<dcmplx>> dD;
        for (const auto& b : D0)
        {
            const size_t n=b.bs->GetNumFunctions();
            chmat_t X(n);
            for (size_t i=0;i<n;i++) {X(i,i)=u(rng); for (size_t j=i+1;j<n;j++) X(i,j)=dcmplx(u(rng),u(rng));}
            dD.push_back({b.irrep, b.bs, X});
        }
        return dD;
    };
    auto apply=[&](std::vector<GroundBlock<dcmplx>> dD, double scale)
    {
        for (auto& b : dD) b.D*=scale;
        auto delta=ChargeDensity::AO_TransitionDensity_Factory<dcmplx>(AsTransition(dD), std::make_shared<Symmetry::Invariant>());
        auto F=K->InducedFock(*delta);
        std::vector<mat_t<dcmplx>> out;
        for (const auto& b : dD) out.push_back(mat_t<dcmplx>(F->Matrix(b.irrep, b.irrep)/scale));
        return out;
    };
    auto relDiff=[](const std::vector<mat_t<dcmplx>>& a, const std::vector<mat_t<dcmplx>>& b)
    {
        double d=0, sc=0;
        for (size_t k=0;k<a.size();k++)
            for (size_t i=0;i<a[k].rows();i++)
                for (size_t j=0;j<a[k].columns();j++) {d=std::max(d, std::abs(a[k](i,j)-b[k](i,j))); sc=std::max(sc, std::abs(b[k](i,j)));}
        return d/sc;
    };
    const auto d1=random(), d2=random();
    const auto K1=apply(d1, 1.0);
    // (a) THE RAW KERNEL.  Through R2 it rode the ground-state collocation, whose D-aware screen has an ABSOLUTE
    // tolerance, so it was scale-dependent (2026-09-29: 2.6e-7 at s = 1e-2, 6.9e-5 at 1e-5, 3.9% at 1e-8) and only
    // printed here.  Since R3 step 3 the transition collocations (B1/B2) take the GEOMETRY-ONLY screen by
    // construction, and the kernel is homogeneous outright (measured 2.5e-15 at every scale): now asserted.
    for (double sc : {1e-2, 1e-5, 1e-8})
    {
        const double r=relDiff(apply(d1, sc), K1);
        std::cout << "[linearity] RAW K[s dD]/s vs K[dD]  s=" << sc << "  rel " << r << std::endl;
        EXPECT_LT(r, 1e-12) << "the raw kernel is not homogeneous at scale " << sc;
    }
    auto d12=d1;
    for (size_t k=0;k<d12.size();k++) d12[k].D+=d2[k].D;
    const auto K12=apply(d12, 1.0), K2=apply(d2, 1.0);
    std::vector<mat_t<dcmplx>> sum;
    for (size_t k=0;k<K1.size();k++) sum.push_back(mat_t<dcmplx>(K1[k]+K2[k]));
    const double add=relDiff(K12, sum);
    std::cout << "[linearity] RAW K[d1+d2] vs K[d1]+K[d2]  rel " << add << std::endl;
    EXPECT_LT(add, 1e-7) << "the kernel is not additive at unit scale";

    // (b) THE OPERATOR GMRES USES (InducedFockMO: rescaled to unit max-norm in, scaled back out) -- exactly
    // homogeneous by construction, at every scale a Krylov vector can have.
    const Response::Reference ref=Response::MakeReference(wf, OccupationConfig{}, {.acrossK=true, .acrossSpin=false},
                                                          std::numeric_limits<double>::quiet_NaN());
    const auto frame=Response::MakeOrbitalFrame(ref, wf);
    const auto rule=std::make_shared<Symmetry::Invariant>();
    Response::BlockPairs x;
    for (size_t b=0;b<ref.NumBlocks();b++)
    {
        const size_t n=ref.NumOrbitals(b);
        Response::cmat_t X(n, n);
        for (size_t i=0;i<n;i++) {X(i,i)=u(rng); for (size_t j=i+1;j<n;j++) {X(i,j)=dcmplx(u(rng),u(rng)); X(j,i)=std::conj(X(i,j));}}
        x.m.push_back(X);
    }
    auto opApply=[&](double sc)
    {
        Response::BlockPairs y=x;
        for (auto& m : y.m) m*=sc;
        Response::BlockPairs F=Response::InducedFockMO(frame, *K, y, rule);
        std::vector<mat_t<dcmplx>> out;
        for (auto& m : F.m) out.push_back(mat_t<dcmplx>(m/sc));
        return out;
    };
    const auto O1=opApply(1.0);
    for (double sc : {1e-2, 1e-5, 1e-8})
    {
        const double r=relDiff(opApply(sc), O1);
        std::cout << "[linearity] OPERATOR (rescaled) s=" << sc << "  rel " << r << std::endl;
        EXPECT_LT(r, 1e-13) << "the response operator is not homogeneous at scale " << sc;
    }
}

//! R2 end to end through the facade: the q = 0 SELF-CONSISTENT response of GPW Si over a Si-p manifold at U = 0
//! (the U_0 case: +U answers zero).  Physics gates that need no oracle: chi is Hermitian, SCREENED
//! (|chi| < |chi0| on the diagonal, same sign -- Hartree + ALDA reduce an insulator's response), and the Krylov
//! solve converged.  UnPol and Pol must agree (the closed shell imposed polarized).
TEST(ResponsePolarizability, GPW_Si_HubbardLinearResponse_Screened)
{
    std::vector<mat_t<dcmplx>> chis;
    for (SpinGroup g : {SpinGroup::UnPolarized, SpinGroup::Polarized})
    {
        auto calc=ConvergedSi(g);
        auto r=calc->HubbardLinearResponse();
        ASSERT_TRUE(r.IsOk()) << r.Error().detail;
        const size_t n=r->labels.size();
        ASSERT_GT(n, 0u);
        for (size_t I=0;I<n;I++)
        {
            EXPECT_LT(r->chi0(I,I).real(), 0.0) << r->labels[I];
            EXPECT_LT(r->chi (I,I).real(), 0.0) << r->labels[I];
            EXPECT_LT(std::abs(r->chi(I,I)), std::abs(r->chi0(I,I))) << r->labels[I] << ": the kernel must SCREEN";
            for (size_t J=0;J<n;J++) EXPECT_NEAR(std::abs(r->chi(I,J)-std::conj(r->chi(J,I))), 0.0, 1e-8) << "chi not Hermitian";
            EXPECT_LE(r->residual[I], 1e-8);
        }
        chis.push_back(r->chi);
    }
    for (size_t I=0;I<chis[0].rows();I++) EXPECT_NEAR(std::abs(chis[0](I,I)-chis[1](I,I)), 0.0, 1e-6*std::abs(chis[0](I,I))) << "Pol != UnPol";
}

//! GATE (b), R2: the self-consistent chi is dn/dalpha -- LR-cDFT's own definition (Cococcioni & de Gironcoli 2005,
//! Timrov §III): converge at alpha = +-a on the Si-p projector (HubbardManifold::alpha, QE's Hubbard_alpha) and
//! difference the manifold occupation.  An INDEPENDENT route: no kernel, no solver, just two SCFs.  The
//! occupation is read through the same projector amplitudes the probe perturbs with (Adjoint/Forward).
namespace {
double ManifoldOccupation(const SolidCalculation& calc)
{
    const auto& wf=ResponseFacadeTests::WF(calc);
    const auto* hub=ResponseFacadeTests::Ham(calc).GetHubbardChannels();
    const Response::Reference ref=Response::MakeReference(wf, OccupationConfig{}, {}, std::numeric_limits<double>::quiet_NaN());
    const auto probe=Response::MakeHubbardProbe(ref, wf, *hub);
    Response::BlockPairs D;                                    // the ground state in its own MO basis: diag(occupation)
    for (const Irrep& ir : wf.GetQNs())
    {
        const auto* os=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(wf.GetOrbitals(ir));
        const size_t n=os->GetNumOrbitals();
        Response::cmat_t X(n, n, dcmplx(0.0));
        size_t i=0;
        for (const auto* o : os->Iterate<Orbitals::TOrbital<dcmplx>>()) {X(i,i)=o->GetOccupation(); i++;}
        D.m.push_back(X);
    }
    return probe.Measure(Symmetry::Invariant(), D)[0].real();
}
} // namespace

TEST(ResponsePolarizability, GPW_Si_Chi_eqFiniteDifferenceCDFT)
{
    const double a=1e-3;                                       // Ha
    auto c0=ConvergedSi(SpinGroup::UnPolarized);
    auto r=c0->HubbardLinearResponse();
    ASSERT_TRUE(r.IsOk()) << r.Error().detail;
    const double chiLR=r->chi(0,0).real();
    auto cp=ConvergedSi(SpinGroup::UnPolarized, +a), cm=ConvergedSi(SpinGroup::UnPolarized, -a);
    const double np=ManifoldOccupation(*cp), nm=ManifoldOccupation(*cm), n0=ManifoldOccupation(*c0);
    const double chiFD=(np-nm)/(2*a);
    std::cout << "[R2 cDFT] n(-a) " << nm << "  n(0) " << n0 << "  n(+a) " << np << "  chi FD " << chiFD
              << "  chi LR " << chiLR << "  rel " << (chiLR-chiFD)/chiFD << std::endl;
    EXPECT_NEAR(chiLR, chiFD, 1e-4*std::fabs(chiFD));
}

//=====================================================================================================
//  §3d STEP 1: +U FROZEN (Q6) and the PERTURBED set (Q10), at U != 0.  Linear-response U is DEFINED with V_Hub
//  held at its ground-state value (Timrov eq 20); the finite-difference LRT cross-check must then hold it too.
//  Si-p at U = 2 eV on site 0 (perturbed by default: it carries U) + site 1's p at U = 0 (a spectator: measured,
//  not perturbed).  Three claims: (1) frozen LR chi == frozen FD chi, every measured row; (2) the freeze MATTERS
//  here -- an UNFROZEN FD chi (fresh +-alpha runs, as R2's gate) differs by far more than the tolerance, so (1)
//  is not vacuous; (3) the run's +U term comes back unfrozen.
//=====================================================================================================
namespace {
//! \a densityMixing: the +-alpha SCFs run on the Kerker/Pulay recipe (the TMO one) instead of SiParams' linear
//! D-mixing, which is not robust on a restart (see below).  \a tol: the FD oracle's own limit (its SCFs converge to
//! Δρ ~ 1e-7).  \a requireRestore: assert the closing unfrozen SCF converged -- on the k211 recipe it sits at the ground
//! state's energy to 10 digits yet never meets ΔE/E < 1e-10 (a criterion floor, not the gate's claim).
void FrozenChiEqFiniteDifference(ivec3_t kmesh, bool densityMixing, double tol, bool requireRestore)
{
    const double U=2.0/27.211386245988, a=1e-3;
    auto c=ConvergedSi(SpinGroup::UnPolarized, 0.0, U, /*spectator*/true, kmesh);
    auto r=c->HubbardLinearResponse();
    ASSERT_TRUE(r.IsOk()) << r.Error().detail;
    ASSERT_EQ(r->labels.size(), 2u);
    EXPECT_EQ(r->perturbed, std::vector<size_t>{0}) << "Q10: only the manifold that carries U is perturbed by default";
    EXPECT_EQ(r->chi.columns(), 1u);
    EXPECT_FALSE(ResponseFacadeTests::Ham(*c).GetHubbardUTarget()->OccupationsFrozen()) << "the freeze was not restored";

    // RELAX 0.2 for the +-alpha SCFs (user's diagnosis, 2026-09-29): SiParams mixes at relax 1.0, which DIIS hides
    // from the seed but not on a restart, where the error vectors are ONE mode and DIIS keeps 2 of them -- the
    // restarted iteration then grew ~1.3x per two steps.  At 0.2 every restart converges (16/22/14) and FD == LR
    // to 1e-6; at 0.4 the unfrozen restore still failed.
    // WHY U = 2 eV, not 4 (measured 2026-09-29): Dudarev's +U ANTI-screens (unfrozen chi -22.5 vs frozen -14.3 at
    // 4 eV), and at 4 eV the UNFROZEN Si SCF sits near that instability -- its convergence depended on the start
    // (the closing restore wandered for 200 iterations).  2 eV keeps the freeze's effect far above the tolerance.
    SCFParams fdp=SiParams(); fdp.NMaxIter=200; fdp.StartingRelaxRo=0.2;
    if (densityMixing) {fdp.PulayDepth=8; fdp.PulayStart=5; fdp.KerkerG0=1.0; fdp.StartingRelaxRo=0.45;}
    auto fd=c->HubbardFiniteDifferenceChi(0, a, fdp);
    ASSERT_TRUE(fd.IsOk()) << fd.Error().details;
    if (requireRestore) EXPECT_TRUE(fd->restored);
    else std::cout << "[step1 LRT] restore " << (fd->restored ? "converged" : "NOT converged (criterion floor; not asserted)") << std::endl;
    EXPECT_FALSE(ResponseFacadeTests::Ham(*c).GetHubbardUTarget()->OccupationsFrozen()) << "the FD run left +U frozen";
    for (size_t I=0;I<2;I++)
    {
        const double lr=r->chi(I,0).real();
        std::cout << "[step1 LRT] " << r->labels[I] << "  chi LR (frozen) " << lr << "  chi FD (frozen) " << fd->chi[I]
                  << "  rel " << (lr-fd->chi[I])/fd->chi[0] << std::endl;
        EXPECT_NEAR(lr, fd->chi[I], tol*std::fabs(fd->chi[0])) << r->labels[I];
    }
    // (2) the unfrozen FD: +U follows the density, so its kernel screens the response differently.
    auto cp=ConvergedSi(SpinGroup::UnPolarized, +a, U, true, kmesh), cm=ConvergedSi(SpinGroup::UnPolarized, -a, U, true, kmesh);
    const double chiUnfrozen=(ManifoldOccupation(*cp)-ManifoldOccupation(*cm))/(2*a);
    std::cout << "[step1 LRT] unfrozen FD chi " << chiUnfrozen << " vs frozen " << fd->chi[0] << std::endl;
    EXPECT_GT(std::fabs(chiUnfrozen-fd->chi[0]), 100*1e-4*std::fabs(fd->chi[0])) << "the freeze made no difference: the gate is vacuous";
}
} // namespace
TEST(ResponsePolarizability, GPW_Si_U2_FrozenChi_eqFiniteDifferenceLRT) {FrozenChiEqFiniteDifference(ivec3_t(1,1,1), false, 1e-5, true);}
//! The same claim on a MULTI-k mesh (w_k = 1/2): the gate that would have caught the missing BZ weight in the
//! transition density (NiO k222, 2026-09-29) -- every other kernel gate is Γ-only, where w = 1.  Measured: LR == FD
//! to 9e-6.  The FD SCFs take the Kerker/Pulay recipe: SiParams' linear D-mixing did not converge the +-alpha runs on
//! this mesh from either start (restart at relax 0.2, or the seed).
TEST(ResponsePolarizability, GPW_Si_k211_U2_FrozenChi_eqFiniteDifferenceLRT) {FrozenChiEqFiniteDifference(ivec3_t(2,1,1), true, 3e-5, false);}


//=====================================================================================================
//  R3 STEP 3: THE KERNEL AT q != 0 (doc/LinearResponsePlan.md §3d).  No finite-difference oracle exists for a
//  (k+q, k) transition density -- it is not a density -- so the claim is structural: the whole kernel (Hartree at
//  G+q + ALDA f_xc) is a HERMITIAN operator on the pairs, <X, K Y> = conj <Y, K X>, because each term is A^† M A
//  with M real (4π/|G+q|^2, w f_xc).  That holds ONLY if every B2 forward/adjoint pair is exact AND the terms
//  route bra and ket the right way round, so a swapped phase, a missing conjugate or a mis-keyed pair breaks it.
//  At a NON-TRIM q (1/3 on a 3x1x1 mesh: complex phases, feedback_complex_type_vs_value), on BOTH XC samplers
//  (the uniform raster's raw pair route, and the Becke point route), polarized as well.  The ground state need not
//  be converged: Hermiticity is an operator identity at any positive ρ0.  (The supercell equivalence, step 4, is
//  what then checks the VALUES.)
//=====================================================================================================
namespace {
std::unique_ptr<SolidCalculation> SiState(ivec3_t kmesh, qcMesh::UnitCellKind xc, SpinGroup g)
{
    FCCUnitCell cell(10.26);
    cell.AddAtom(14, {0,0,0});
    cell.AddAtom(14, {0.25,0.25,0.25});
    Lattice_3D lat(cell, kmesh);
    auto mol=std::shared_ptr<const BasisSet::Real_BS>(
        BasisSet::Gaussian::Factory(BasisSet::Gaussian::BasisSetData::SIPP_SR, &cell,
                                    BasisSet::Gaussian::Engine::MnD, BasisSet::Gaussian::Angular::Cartesian));
    SCFParams par=SiParams();
    par.NMaxIter=12;                                         // not converged, and need not be (see above)
    SolidCalcOptions o{.Nelec=8, .multiplicity = g==SpinGroup::Polarized ? 1 : 0, .species={{"Si",4}},
                       .densityEcut=20.0, .forceComplex=true};
    o.xcMesh.cellKind=xc;
    return std::make_unique<SolidCalculation>(lat, mol, o, par);
}

void KernelIsHermitianAtNonTrimQ(qcMesh::UnitCellKind xc, SpinGroup g)
{
    auto calc=SiState(ivec3_t(3,1,1), xc, g);
    const auto& wf=ResponseFacadeTests::WF(*calc);
    auto& H=ResponseFacadeTests::Ham(*calc);
    const auto D0cd=wf.GetChargeDensity();
    auto K=H.MakeResponseKernel(&calc->Basis(), D0cd.get());
    const Response::Reference ref=Response::MakeReference(wf, OccupationConfig{}, {.acrossK=true, .acrossSpin=false},
                                                          std::numeric_limits<double>::quiet_NaN());
    auto qs=ref.QMesh(ivec3_t(3,1,1));
    ASSERT_TRUE(qs.IsOk()) << qs.Error().detail;
    std::shared_ptr<const Symmetry::Lattice_3D::MeshShift> rule;
    for (const auto& q : qs.Value()) if (q.Steps().x==1) rule=std::make_shared<const Symmetry::Lattice_3D::MeshShift>(q);
    ASSERT_TRUE(rule);
    // The (k+q, k) pairs over the wave function's blocks, and two random δD sets on them.
    struct Blk {Irrep ir; const Hamiltonian::cobs_t* bs;};
    std::vector<Blk> blocks;
    for (const Irrep& ir : wf.GetQNs())
    {
        const auto* os=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(wf.GetOrbitals(ir));
        blocks.push_back({ir, dynamic_cast<const Hamiltonian::cobs_t*>(os->GetBasisSet())});
    }
    std::mt19937 rng(20260929);
    std::uniform_real_distribution<double> u(-0.05, 0.05);
    auto random=[&]()
    {
        std::vector<TransitionBlock<dcmplx>> out;
        for (const auto& ket : blocks)
            for (const auto& bra : blocks)
            {
                if (bra.ir.ms!=ket.ir.ms || !rule->Couples(*bra.ir.sym, *ket.ir.sym)) continue;
                mat_t<dcmplx> X(bra.bs->GetNumFunctions(), ket.bs->GetNumFunctions());
                for (size_t i=0;i<X.rows();i++) for (size_t j=0;j<X.columns();j++) X(i,j)=dcmplx(u(rng),u(rng));
                out.push_back({bra.ir, ket.ir, bra.bs, ket.bs, X});
            }
        return out;
    };
    const auto X=random(), Y=random();
    ASSERT_EQ(X.size(), blocks.size()) << "every ket block needs exactly one (k+q) partner";
    auto FX=K->InducedFock(*ChargeDensity::AO_TransitionDensity_Factory<dcmplx>(X, rule));
    auto FY=K->InducedFock(*ChargeDensity::AO_TransitionDensity_Factory<dcmplx>(Y, rule));
    auto pair=[](const std::vector<TransitionBlock<dcmplx>>& A, const Hamiltonian::TransitionFock<dcmplx>& F)
    {
        dcmplx s=0.0;
        for (const auto& p : A)
        {
            const mat_t<dcmplx> f=F.Matrix(p.bra, p.ket);
            for (size_t i=0;i<f.rows();i++) for (size_t j=0;j<f.columns();j++) s+=std::conj(p.dD(i,j))*f(i,j);
        }
        return s;
    };
    const dcmplx xy=pair(X, *FY), yx=pair(Y, *FX), xx=pair(X, *FX);
    std::cout << "[R3 kernel q=1/3] <X,KY> " << xy << "  conj<Y,KX> " << std::conj(yx) << "  rel "
              << std::abs(xy-std::conj(yx))/std::abs(xy) << "   <X,KX> " << xx << std::endl;
    EXPECT_GT(std::abs(xy), 1e-6) << "a vacuous comparison";
    EXPECT_LT(std::abs(xy-std::conj(yx)), 1e-11*std::abs(xy)) << "the q != 0 kernel is not Hermitian";
    EXPECT_LT(std::fabs(xx.imag()), 1e-11*std::abs(xx)) << "<X,KX> must be real";
}
} // namespace

TEST(ResponseKernel, GPW_Si_k311_q13_KernelIsHermitian_Raster_UnPol) {KernelIsHermitianAtNonTrimQ(qcMesh::UnitCellKind::Uniform, SpinGroup::UnPolarized);}
TEST(ResponseKernel, GPW_Si_k311_q13_KernelIsHermitian_Raster_Pol)   {KernelIsHermitianAtNonTrimQ(qcMesh::UnitCellKind::Uniform, SpinGroup::Polarized);}
TEST(ResponseKernel, GPW_Si_k311_q13_KernelIsHermitian_Becke_UnPol)  {KernelIsHermitianAtNonTrimQ(qcMesh::UnitCellKind::Becke,   SpinGroup::UnPolarized);}

//! §3d FINDING 5, CLOSED (R3 step 3): an IMPOSED Γ run on the uniform raster -- the T3 stream fold and the raster
//! star-average both armed -- gives the FREE run's chi.  Before step 3 the transition density rode the ground-state
//! composite, whose leaves folded δD and star-averaged δρ under a group the one-site perturbation breaks, so
//! HubbardLinearResponse REFUSED this configuration.  The pair route has no fold in it, so the refusal is gone.
TEST(ResponsePolarizability, GPW_Si_ImposedGamma_Raster_eqFreeChi)
{
    auto f=ConvergedSi(SpinGroup::UnPolarized);
    auto i=ConvergedSi(SpinGroup::UnPolarized, 0.0, 0.0, false, ivec3_t(1,1,1), /*impose*/true);
    auto rf=f->HubbardLinearResponse(), ri=i->HubbardLinearResponse();
    ASSERT_TRUE(rf.IsOk()) << rf.Error().detail;
    ASSERT_TRUE(ri.IsOk()) << ri.Error().detail;
    const double cf=rf->chi(0,0).real(), ci=ri->chi(0,0).real(), c0f=rf->chi0(0,0).real(), c0i=ri->chi0(0,0).real();
    std::cout << "[finding 5] chi0 free " << c0f << " imposed " << c0i << "   chi free " << cf << " imposed " << ci
              << "  rel " << (ci-cf)/cf << std::endl;
    // MEASURED 2026-09-29: chi0 5e-6, chi 9.5e-6 relative.  That is the two GROUND STATES' agreement (E differs by
    // 4e-7 Ha; chi0 has no kernel in it and already differs by 5e-6), not the kernel.  A symmetrized δρ -- the
    // one-site probe averaged with its image on the other Si -- changes the kernel's INPUT at O(1), so 5e-5 keeps
    // the gate far below the defect it guards against while clear of the ground-state floor.
    EXPECT_NEAR(c0i, c0f, 5e-5*std::fabs(c0f));
    EXPECT_NEAR(ci,  cf,  5e-5*std::fabs(cf)) << "an imposed run's response differs from the free run's";
}

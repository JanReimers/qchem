// File: Calculation/tests/ForceFacade.C  The facade-side pieces of the analytic force (doc/GeometryOptimizationPlan.md §4-§6, §13),
// with no Hamiltonian-hub change: the energy-weighted density W (the Pulay term's input), and the PER-TERM FINITE DIFFERENCE AT FIXED D
// -- one term's energy at R ± h with the AO density frozen, no SCF -- gated on the nuclear attraction, whose analytic gradient libcint
// gives (int1e_ipnuc + int1e_iprinv, PG_LibCint::NuclearGradient), and on the kinetic term (int1e_ipkin).
#include <gtest/gtest.h>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <memory>
#include <vector>
#include "forward.H"
import qchem.MoleculeCalculation;
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.BasisSet;
import qchem.Hamiltonian;
import qchem.Symmetry.Spin;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;
import qchem.BasisSet.AoShellSource;   // the AD kernel's shells
import qchem.Mesh.Product;             // the AD kernel's one-centre mesh template
import qchem.AD.XCKernel;              // the Enzyme AD spike's kernel

using namespace qchem;
using StructureData::GetMolecule;
using LibCintEval = BasisSet::Gaussian::Evaluators::PG_LibCint::NR_Evaluator;
using PGData      = BasisSet::Gaussian::Evaluators::PG_Cart_MnD::PGData;

//! The friend (forward.H): reaches the facade's basis, Hamiltonian and carry -- what a TEST of one term needs, never a client.
class ForceFacadeTests
{
public:
    //! Σ over the blocks of Tr(W S): the orthonormality identity W must satisfy, Σ_levels f ε.
    static double TraceWS(const MoleculeCalculation& c)
    {
        const auto W=c.MakeEnergyWeightedDensity();
        size_t b=0; double t=0.0;
        for (Spin s : SpinIrreps(c.itsSetup->opts.spin)) { (void)s;
            for (auto ibs : c.itsBasis->Iterate<BasisSet::Real_OIBS>())
            {
                const auto& S=ibs->Overlap();
                for (size_t i=0;i<S.rows();++i) for (size_t j=0;j<S.rows();++j) t+=W[b](i,j)*S(i,j);
                ++b;
            } }
        EXPECT_EQ(b, W.size());
        return t;
    }
    static MoleculeCalculation::Carry Carry(const MoleculeCalculation& c) {return c.MakeCarry();}
    //! THE FIXED-D HELPER: every term's energy at structure \a at, with the carried AO density \a D frozen on that point's basis.
    static std::vector<Hamiltonian::TermEnergy> TermEnergiesAt(std::shared_ptr<const MoleculeSetup> setup, const Structure& at,
                                                               const MoleculeCalculation::Carry& D)
    {
        MoleculeCalculation c(MoleculeCalculation::BuildOnly{}, setup, at, nullptr, nullptr);
        c.BuildGraph();
        const auto cd=c.MakeCarriedDensity(D);
        c.itsHam->RefreshForDensity(c.itsBasis, cd.get());
        return c.itsHam->MakeTermEnergies(cd.get());
    }
    //! Σ_ij D_ij dM_ij/dR_A by libcint on the converged point's (single, Cartesian) block: \a which 0 = kinetic, 1 = nuclear.
    static std::vector<rvec3_t> Contracted(const MoleculeCalculation& c, const MoleculeCalculation::Carry& D, int which)
    {
        size_t b=0; std::vector<rvec3_t> g;
        for (auto ibs : c.itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            const auto* pg=dynamic_cast<const PGData*>(ibs);   // a test may cheat: the Cartesian block IS its PGData
            EXPECT_NE(pg, nullptr);
            LibCintEval lc(*pg, c.itsStructure.get());
            const auto dM = which==0 ? lc.KineticGradient() : lc.NuclearGradient();
            if (g.empty()) g.assign(dM.size(), rvec3_t(0,0,0));
            for (size_t A=0;A<dM.size();++A)
            {
                double x[3]={0,0,0};
                for (int k=0;k<3;++k) for (size_t i=0;i<dM[A][k].rows();++i) for (size_t j=0;j<dM[A][k].columns();++j) x[k]+=D.D[b](i,j)*dM[A][k](i,j);
                g[A]+=rvec3_t(x[0],x[1],x[2]);
            }
            ++b;
        }
        return g;
    }
    static const Structure& Structure_(const MoleculeCalculation& c) {return *c.itsStructure;}
    //! THE HAMILTONIAN'S FORCE (GeometryOptimizationPlan §13) at the carried density \a D with the energy-weighted density \a W.
    static rvec3vec_t Forces(const MoleculeCalculation& c, const MoleculeCalculation::Carry& D, const MoleculeCalculation::Carry& W)
    {
        const auto cd=c.MakeCarriedDensity(D), wd=c.MakeCarriedDensity(W);
        return c.itsHam->At(c.itsBasis, cd.get()).Forces(wd.get(), c.itsStructure.get());
    }
    static std::vector<std::string> MissingGradients(const MoleculeCalculation& c) {return c.itsHam->MissingGradients(c.itsBasis, c.itsStructure.get());}
    static size_t NumBlocks(const MoleculeCalculation& c) {return c.itsBasis->GetNumIBS();}
    //! EVERY TERM'S FORCE at structure \a at with \a D frozen (MakeTermEnergies' twin, index for index): TermEnergiesAt's analytic partner.
    static std::vector<Hamiltonian::TermForce> TermForcesAt(std::shared_ptr<const MoleculeSetup> setup, const Structure& at,
                                                            const MoleculeCalculation::Carry& D)
    {
        MoleculeCalculation c(MoleculeCalculation::BuildOnly{}, setup, at, nullptr, nullptr);
        c.BuildGraph();
        const auto cd=c.MakeCarriedDensity(D);
        c.itsHam->RefreshForDensity(c.itsBasis, cd.get());
        return c.itsHam->MakeTermForces(c.itsBasis, cd.get(), c.itsStructure.get());
    }
    static MoleculeCalculation::Carry EnergyWeighted(const MoleculeCalculation& c) {return {c.MakeEnergyWeightedDensity()};}
    //! The whole energy breakdown at structure \a at with \a D frozen: for terms MakeTermEnergies leaves out (the whole-system HF ones).
    static EnergyBreakdown TotalEnergyAt(std::shared_ptr<const MoleculeSetup> setup, const Structure& at, const MoleculeCalculation::Carry& D)
    {
        MoleculeCalculation c(MoleculeCalculation::BuildOnly{}, setup, at, nullptr, nullptr);
        c.BuildGraph();
        const auto cd=c.MakeCarriedDensity(D);
        c.itsHam->RefreshForDensity(c.itsBasis, cd.get());
        return c.itsHam->GetTotalEnergy(cd.get());
    }
    //! THE PULAY SUM, Σ_blocks Tr(W ∂S/∂R) by each block's overlap pullback, in Forces' sign (added to F): for a model whose Forces()
    //! refuses on a term that is identically zero (Xα at α = 0).
    static rvec3vec_t Pulay(const MoleculeCalculation& c, const MoleculeCalculation::Carry& W)
    {
        rvec3vec_t F(c.itsStructure->GetNumAtoms(), rvec3_t(0,0,0));
        size_t b=0;
        for (Spin s : SpinIrreps(c.itsSetup->opts.spin)) { (void)s;
            for (auto ibs : c.itsBasis->Iterate<BasisSet::Real_OIBS>())
            {
                const auto g=ibs->MakeOverlapPullback(W.D[b++], c.itsStructure.get());
                EXPECT_TRUE(bool(g));
                if (g) for (size_t A=0;A<F.size();A++) F[A]+=(*g)[A];
            } }
        EXPECT_EQ(b, W.D.size());
        return F;
    }
    static MoleculeCalculation::Carry Zero(const MoleculeCalculation::Carry& like)
    {
        MoleculeCalculation::Carry z{like.D};
        for (auto& b : z.D) b=hmat_t<double>(b.rows());   // Blaze zero-initialises
        return z;
    }
    //! The term's H at the point \a c with the carried \a D frozen: the sum of the blocks whose terms report \a entries (one block: the
    //! unblocked molecular basis).
    static hmat_t<double> TermMatrix(const MoleculeCalculation& c, const MoleculeCalculation::Carry& D, const std::string& entries)
    {
        MoleculeCalculation b(MoleculeCalculation::BuildOnly{}, c.itsSetup, *c.itsStructure, nullptr, nullptr);
        b.BuildGraph();
        const auto cd=b.MakeCarriedDensity(D);
        b.itsHam->RefreshForDensity(b.itsBasis, cd.get());
        const auto es=b.itsHam->MakeTermEnergies(cd.get());
        const Spin s=SpinIrreps(b.itsSetup->opts.spin).front();
        for (auto ibs : b.itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            const auto Hs=b.itsHam->MakeTermMatrices(ibs, s, cd.get());
            hmat_t<double> H(ibs->GetNumFunctions());   // Blaze zero-initialises
            for (size_t t=0;t<es.size();t++) if (es[t].entries==entries) H+=Hs[t];
            return H;   // one block (checked by the caller against the kernel's nbf)
        }
        return {};
    }
    //! The AD kernel's data copied out of the converged point \a c: its positions' order, the (sized) mesh's one-centre template, the
    //! Cartesian shells of its single block, and the carried total D.
    static AD::XCKernelData XCData(const MoleculeCalculation& c, const MoleculeCalculation::Carry& D)
    {
        AD::XCKernelData d;
        std::vector<rvec3_t> R; c.itsStructure->ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
        d.natom=R.size();
        const qcMesh::MeshParams& mp=c.itsSetup->mesh;
        d.beckeOrder=mp.beckeOrder;
        const qcMesh::Mesh t=qcMesh::ProductMesh(qcMesh::MakeRadial(mp), qcMesh::MakeAngular(mp));
        for (size_t q=0;q<t.size();q++) { d.u.push_back(t.Points()[q].x); d.u.push_back(t.Points()[q].y); d.u.push_back(t.Points()[q].z); d.w.push_back(t.Weights()[q]); }
        size_t nblock=0;
        for (auto ibs : c.itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            ++nblock;
            d.nbf=ibs->GetNumFunctions();
            const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(ibs);
            EXPECT_NE(src, nullptr);
            for (const auto& sh : src->GetAoShells())
            {
                AD::XCShell x;
                double best=1e300;
                for (size_t A=0;A<R.size();A++) if (norm(R[A]-sh.center)<best) { best=norm(R[A]-sh.center); x.atom=A; }
                EXPECT_LT(best, 1e-12) << "a shell sits on an atom";
                x.offset=sh.offset;
                const auto mono=sh.rep->Monomials();
                EXPECT_EQ(mono.size(), sh.nComponents()) << "a Cartesian shell (the kernel's scope)";
                for (const auto& m : mono) x.powers.push_back({m.n, m.l, m.m});
                for (size_t k=0;k<sh.norm.size();k++) x.norms.push_back(sh.norm[k]);
                for (size_t k=0;k<sh.exponents.size();k++) { x.exponents.push_back(sh.exponents[k]); x.coefficients.push_back(sh.coefficients[k]); }
                d.shells.push_back(std::move(x));
            }
        }
        EXPECT_EQ(nblock, 1u) << "the kernel's scope: one unblocked block";
        EXPECT_EQ(D.D.size(), 1u);
        for (size_t i=0;i<d.nbf;i++) for (size_t j=0;j<d.nbf;j++) d.D.push_back(D.D[0](i,j));
        return d;
    }
    static std::vector<double> AOValues(const MoleculeCalculation& c, const rvec3_t& r)
    {
        for (auto ibs : c.itsBasis->Iterate<BasisSet::Real_OIBS>())
        {
            const auto v=(*ibs)(r);
            std::vector<double> out(v.size());
            for (size_t i=0;i<v.size();i++) out[i]=v[i];
            return out;
        }
        return {};
    }
};

namespace
{
SCFParams Tight() {SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-9; return p;}
double Entry(const std::vector<Hamiltonian::TermEnergy>& es, const std::string& name)
{
    double e=0.0; bool found=false;
    for (const auto& t : es) if (t.entries==name) { e+=t.E; found=true; }
    EXPECT_TRUE(found) << "no term reports " << name;
    return e;
}
double C(const rvec3_t& v, int k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
}

TEST(ForceFacade, TheEnergyWeightedDensityTracesToTheOccupiedEigenvalues)
{
    for (auto model : {Model::HF, Model::LDA})
    {
        MoleculeCalculation c(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=model}), Tight());
        ASSERT_TRUE(c.IsConverged());
        double sum=0.0; for (const auto& l : c.GetLevels()) sum+=l.occupation*l.e;
        EXPECT_NEAR(ForceFacadeTests::TraceWS(c), sum, 1e-9*std::fabs(sum)) << "Tr(W S) = Σ f ε";
    }
}

TEST(ForceFacade, ATermsEnergyAtFixedDMovesAsItsAnalyticGradientSays)
{
    // All-electron LDA water, unblocked Cartesian DZVP: Ven and the kinetic term at R ± h with the converged D frozen, against
    // Σ D dV/dR and Σ D dK/dR from libcint.  The finite difference sees the basis centres AND (for Ven) the nuclei move.
    const Molecule h2o=GetMolecule("H2O");
    const auto setup=MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    // The kinetic TERM's energy is ½ Tr(D <-∇²>) (the project's convention: the matrices carry <-∇²>, the ½ is the Hamiltonian's), so its
    // gradient is ½ of the evaluator's: measured ratio 2.000000 before the ½ was applied (2026-11-07).
    auto gK=ForceFacadeTests::Contracted(c, D, 0); for (auto& g : gK) g=0.5*g;
    const auto gV=ForceFacadeTests::Contracted(c, D, 1);
    std::vector<rvec3_t> R0; h2o.ForEachSite([&](int, const rvec3_t& r, bool){ R0.push_back(r); });
    ASSERT_EQ(gV.size(), R0.size());
    const double h=1e-4;
    rvec3_t sV(0,0,0), sK(0,0,0);
    for (size_t A=0;A<R0.size();++A)
    {
        sV+=gV[A]; sK+=gK[A];
        for (int x=0;x<3;++x)
        {
            auto at=[&](double s){ auto R=R0; (x==0 ? R[A].x : x==1 ? R[A].y : R[A].z)+=s; std::shared_ptr<Structure> st=h2o.Clone(); st->MoveSites(R);
                                    return ForceFacadeTests::TermEnergiesAt(setup, *st, D); };
            const auto p=at(h), m=at(-h);
            const double fV=(Entry(p,"Een")-Entry(m,"Een"))/(2*h), fK=(Entry(p,"Kinetic")-Entry(m,"Kinetic"))/(2*h);
            EXPECT_NEAR(C(gV[A],x), fV, 2e-6*std::max(1.0, std::fabs(fV))) << "dEen/dR atom " << A << " axis " << x;
            EXPECT_NEAR(C(gK[A],x), fK, 2e-6*std::max(1.0, std::fabs(fK))) << "dEkin/dR atom " << A << " axis " << x;
        }
    }
    EXPECT_LT(norm(sK), 1e-8) << "the kinetic term is translation invariant";
    EXPECT_GT(norm(gV[0]), 1e-2) << "the oxygen's Een gradient is not trivially zero";
}

// ============================================================================================================================
// THE ENZYME AD SPIKE (doc/GeometryOptimizationPlan.md, "Automatic differentiation spike"; lane B, 2026-11-08).  The Becke-mesh delta
// XC energy at fixed D on water, as qchem.AD.XCKernel's pure kernel of (R, D): the kernel IS the term (its energy and basis values
// equal the library's), its AD gradient in R IS the per-term finite difference (the library term's energy at R +- h with D frozen),
// it sums to zero over the atoms, and its AD gradient in D IS the term's H (pin 32).
// ============================================================================================================================
namespace
{
std::vector<double> Positions(const Structure& st) {std::vector<double> x; st.ForEachSite([&](int, const rvec3_t& r, bool){ x.push_back(r.x); x.push_back(r.y); x.push_back(r.z); }); return x;}
//! Peak RSS since the last reset, and the reset (Linux: /proc/self/clear_refs "5" restarts VmHWM at the current RSS).
long StatusKB(const std::string& key)
{
    std::ifstream f("/proc/self/status"); std::string k; long v=0, out=-1; std::string unit;
    while (f>>k) { if (k==key+":") { f>>v; out=v; break; } std::getline(f,unit); }
    return out;
}
void ResetPeak() {std::ofstream("/proc/self/clear_refs")<<"5";}
}

TEST(ADSpike, TheKernelIsTheDeltaXCTermAtFixedD)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const AD::XCKernelData d=ForceFacadeTests::XCData(c, D);
    const std::vector<double> R=Positions(ForceFacadeTests::Structure_(c));
    for (const rvec3_t& r : {rvec3_t(0.1,0.2,0.3), rvec3_t(-0.8,1.1,-0.4), rvec3_t(1.5,-0.3,0.9), rvec3_t(0,1.43,-0.88)})
    {
        const auto lib=ForceFacadeTests::AOValues(c, r), ker=AD::XCKernelAOValues(d, R, r);
        ASSERT_EQ(lib.size(), ker.size());
        for (size_t i=0;i<lib.size();i++) EXPECT_NEAR(ker[i], lib[i], 1e-13*std::max(1.0,std::fabs(lib[i]))) << "basis function " << i;
    }
    const double Elib=Entry(ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D), "Exc");
    const double Eker=AD::XCEnergy(d, R);
    std::cout<<"[AD spike] E_xc kernel "<<std::setprecision(14)<<Eker<<"  term "<<Elib<<"  diff "<<Eker-Elib<<std::endl;
    EXPECT_NEAR(Eker, Elib, 1e-11*std::fabs(Elib)) << "the kernel is the production delta XC term";
}

TEST(ADSpike, ItsADForceIsThePerTermFiniteDifferenceAndSumsToZero)
{
    const Molecule h2o=GetMolecule("H2O");
    const auto setup=MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const AD::XCKernelData d=ForceFacadeTests::XCData(c, D);
    const std::vector<double> R=Positions(ForceFacadeTests::Structure_(c));
    const auto g=AD::XCEnergyGradientR(d, R);
    // COST: one energy, one gradient (R), one gradient (D), wall time on this process's OpenMP threads; TAPE: the gradient's peak RSS
    // above the RSS it started at.
    using clk=std::chrono::steady_clock;
    auto t0=clk::now(); const double E=AD::XCEnergy(d, R); auto t1=clk::now();
    ResetPeak(); const long rss0=StatusKB("VmRSS");
    auto t2=clk::now(); const auto g2=AD::XCEnergyGradientR(d, R); auto t3=clk::now();
    const long tapeR=StatusKB("VmHWM")-rss0;
    ResetPeak(); const long rss1=StatusKB("VmRSS");
    auto t4=clk::now(); const auto gD=AD::XCEnergyGradientD(d, R); auto t5=clk::now();
    const long tapeD=StatusKB("VmHWM")-rss1;
    const double tE=std::chrono::duration<double>(t1-t0).count(), tR=std::chrono::duration<double>(t3-t2).count(), tD=std::chrono::duration<double>(t5-t4).count();
    std::cout<<"[AD spike] "<<d.w.size()*d.natom<<" mesh points, "<<d.nbf<<" basis functions; E "<<E<<"\n[AD spike] energy "<<tE<<" s, dE/dR "<<tR<<" s ("
             <<tR/tE<<" x), dE/dD "<<tD<<" s ("<<tD/tE<<" x); tape (peak RSS above start) dE/dR "<<tapeR/1024.0<<" MB, dE/dD "<<tapeD/1024.0<<" MB"<<std::endl;
    for (size_t i=0;i<R.size();i++) EXPECT_EQ(g[i], g2[i]) << "the gradient is deterministic (block sums in block order)";
    // THE PER-TERM FINITE DIFFERENCE: the library's delta XC term at R +- h with D frozen (basis centres, mesh points and partition move).
    const double h=1e-4;
    rvec3_t sum(0,0,0);
    for (size_t A=0;A<d.natom;++A)
    {
        sum+=rvec3_t(g[3*A], g[3*A+1], g[3*A+2]);
        for (int x=0;x<3;++x)
        {
            auto at=[&](double s){ std::vector<rvec3_t> Rv; for (size_t B=0;B<d.natom;B++) Rv.push_back(rvec3_t(R[3*B],R[3*B+1],R[3*B+2]));
                                    (x==0 ? Rv[A].x : x==1 ? Rv[A].y : Rv[A].z)+=s; std::shared_ptr<Structure> st=h2o.Clone(); st->MoveSites(Rv);
                                    return Entry(ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *st, D), "Exc"); };
            const double fd=(at(h)-at(-h))/(2*h);
            std::cout<<"[AD spike] dExc/dR atom "<<A<<" axis "<<x<<": AD "<<std::setprecision(10)<<g[3*A+x]<<"  FD "<<fd<<"  diff "<<g[3*A+x]-fd<<std::endl;
            EXPECT_NEAR(g[3*A+x], fd, 2e-6*std::max(1.0, std::fabs(fd))) << "dExc/dR atom " << A << " axis " << x;
        }
    }
    std::cout<<"[AD spike] Σ_A dExc/dR_A = "<<sum<<std::endl;
    EXPECT_LT(norm(sum), 1e-10) << "translation invariance holds to rounding: the grid response (moving points and partition) is in";
}

TEST(ADSpike, ItsDGradientIsTheTermsH)
{
    const auto setup=MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=Model::LDA});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const AD::XCKernelData d=ForceFacadeTests::XCData(c, D);
    const auto g=AD::XCEnergyGradientD(d, Positions(ForceFacadeTests::Structure_(c)));
    const hmat_t<double> H=ForceFacadeTests::TermMatrix(c, D, "Exc");
    ASSERT_EQ(H.rows(), d.nbf);
    double worst=0.0, largest=0.0;
    for (size_t i=0;i<d.nbf;i++) for (size_t j=0;j<d.nbf;j++)
    {
        worst=std::max(worst, std::fabs(g(i,j)-H(i,j))); largest=std::max(largest, std::fabs(H(i,j)));
        EXPECT_NEAR(g(i,j), H(i,j), 1e-11) << "dE/dD_" << i << j << " against the term's H";
    }
    std::cout<<"[AD spike] pin 32: max|dE/dD - H| = "<<worst<<" (max|H| "<<largest<<")"<<std::endl;
}

// ============================================================================================================================
// THE HAMILTONIAN'S FORCE (GeometryOptimizationPlan §13; lane B 2026-11-08), layer by layer on the one-electron model E1 (Kinetic + Ven +
// IonIon, no electron-electron term), where every term has its force: (1) at D = 0 and W = 0 it is the ion-ion force alone; (2) at the
// converged D with W = 0 it is the terms' Hellmann-Feynman sum, the per-term finite difference at fixed D; (3) with W it is the whole force,
// the finite difference of the RE-CONVERGED energy (the Pulay term is what makes the fixed-D sum into that).  MissingGradients() is empty
// for E1 and names the electron-electron terms for HF and LDA.
// ============================================================================================================================
namespace
{
std::shared_ptr<Structure> MovedH2O(const Molecule& h2o, size_t A, int x, double s)
{
    std::vector<rvec3_t> R; h2o.ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
    (x==0 ? R[A].x : x==1 ? R[A].y : R[A].z)+=s;
    std::shared_ptr<Structure> st=h2o.Clone(); st->MoveSites(R);
    return st;
}
}

TEST(HamiltonianForces, OnBareNucleiTheForceIsTheIonIonForce)
{
    const Molecule h2o=GetMolecule("H2O");
    const auto setup=MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::E1});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    EXPECT_TRUE(ForceFacadeTests::MissingGradients(c).empty()) << "E1 has every term's force";
    const auto D0=ForceFacadeTests::Zero(ForceFacadeTests::Carry(c));
    const rvec3vec_t F=ForceFacadeTests::Forces(c, D0, D0);   // D = 0, W = 0: the electrons contribute nothing
    // E_nn is a closed form, so the FD is limited by its own h²/6 |E‴| ~ 2e-10 at h = 1e-5: 1e-8.
    const double h=1e-5;
    for (size_t A=0;A<F.size();A++) for (int x=0;x<3;x++)
    {
        const double fd=(Entry(ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *MovedH2O(h2o,A,x,h), D0), "Enn")
                        -Entry(ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *MovedH2O(h2o,A,x,-h), D0), "Enn"))/(2*h);
        EXPECT_NEAR(C(F[A],x), -fd, 1e-8) << "atom " << A << " axis " << x;
    }
}

TEST(HamiltonianForces, AtFixedDTheTermsForceIsThePerTermFiniteDifference)
{
    const Molecule h2o=GetMolecule("H2O");
    const auto setup=MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::E1});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const rvec3vec_t F=ForceFacadeTests::Forces(c, D, ForceFacadeTests::Zero(D));   // W = 0: the terms' sum alone
    // The fixed-D energy is smooth in R; at h = 1e-4 the central difference's h²/6 |E‴| is ≲ 1e-8, so 2e-6 (relative above 1).
    const double h=1e-4;
    rvec3_t sum(0,0,0);
    for (size_t A=0;A<F.size();A++)
    {
        sum+=F[A];
        for (int x=0;x<3;x++)
        {
            auto E=[&](double s){ const auto es=ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *MovedH2O(h2o,A,x,s), D);
                                  return Entry(es,"Kinetic")+Entry(es,"Een")+Entry(es,"Enn"); };
            const double fd=(E(h)-E(-h))/(2*h);
            EXPECT_NEAR(C(F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd))) << "atom " << A << " axis " << x;
        }
    }
    EXPECT_LT(norm(sum), 1e-8) << "translation invariance (the basis, the nuclei and the ions all move together)";
}

TEST(HamiltonianForces, WithThePulayTermItIsTheFiniteDifferenceOfTheReconvergedEnergy)
{
    const Molecule h2o=GetMolecule("H2O");
    const auto setup=MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::E1});
    MoleculeCalculation c(setup, Tight());
    ASSERT_TRUE(c.IsConverged());
    const rvec3vec_t F=ForceFacadeTests::Forces(c, ForceFacadeTests::Carry(c), ForceFacadeTests::EnergyWeighted(c));
    // Each displaced point is a full SCF (MinΔρ/N 1e-9: its energy error is second order, ≲ 1e-12).  h = 1e-4 puts the FD's h²/6 |E‴| at
    // ≲ 1e-8, so 2e-6 (relative above 1) leaves room for both.
    const double h=1e-4;
    for (size_t A=0;A<F.size();A++) for (int x=0;x<3;x++)
    {
        auto E=[&](double s){ MoleculeCalculation m(MoleculeSetup::Resolve(dynamic_cast<const Molecule&>(*MovedH2O(h2o,A,x,s)), {.basis="dzvp", .model=Model::E1}), Tight());
                              EXPECT_TRUE(m.IsConverged()); return m.Energy(); };
        const double fd=(E(h)-E(-h))/(2*h);
        EXPECT_NEAR(C(F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd))) << "atom " << A << " axis " << x;
    }
}

TEST(HamiltonianForces, MissingGradientsNamesTheTermsWithoutAForce)
{
    // LDA: the XC term has no force yet.  (HF has every term's force since the HF Coulomb and exchange landed: see its gates below.)
    for (auto model : {Model::LDA})
    {
        MoleculeCalculation c(MoleculeSetup::Resolve(GetMolecule("H2O"), {.basis="dzvp", .model=model}), Tight());
        ASSERT_TRUE(c.IsConverged());
        const auto missing=ForceFacadeTests::MissingGradients(c);
        EXPECT_FALSE(missing.empty()) << "the XC term has no force yet";
        for (const auto& m : missing) { EXPECT_NE(m.find("no geometry gradient"), std::string::npos) << m; std::cout<<"[forces] missing: "<<m<<std::endl; }
        const auto D=ForceFacadeTests::Carry(c);
        EXPECT_THROW((void)ForceFacadeTests::Forces(c, D, D), std::logic_error) << "Forces refuses what MissingGradients names";
    }
}

// ============================================================================================================================
// THE PER-TERM FORCE (MakeTermForces, MakeTermEnergies' derivative twin) and the DUNLAP HARTREE term's force (§13.2; lane B 2026-11-08).
// ============================================================================================================================
TEST(HamiltonianForces, ThePerTermForcesSumToTheForceWithoutItsPulayTerm)
{
    // E1 water: every term has its force, so Σ_t F_t is Forces at W = 0, term order is MakeTermEnergies', and nothing is missing.
    const Molecule h2o=GetMolecule("H2O");
    MoleculeCalculation c(MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::E1}), Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const auto Ft=ForceFacadeTests::TermForcesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    const auto Et=ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    ASSERT_EQ(Ft.size(), Et.size()) << "one force per term energy";
    const rvec3vec_t F=ForceFacadeTests::Forces(c, D, ForceFacadeTests::Zero(D));
    rvec3vec_t sum(F.size(), rvec3_t(0,0,0));
    for (const auto& t : Ft) { ASSERT_TRUE(t.missing.empty()) << t.missing; ASSERT_EQ(t.F.size(), F.size()); for (size_t A=0;A<F.size();A++) sum[A]+=t.F[A]; }
    for (size_t A=0;A<F.size();A++) EXPECT_LT(norm(sum[A]-F[A]), 1e-12) << "atom " << A;
}

TEST(HamiltonianForces, WithTheDunlapHartreeTheForceIsTheFiniteDifferenceOfTheReconvergedEnergy_Long)
{
    // Xα at α = 0: the XC term adds EXACTLY nothing (energy and potential ∝ α), so the Hamiltonian is Kinetic + Ven + IonIon + the default
    // (Dunlap) Hartree, and the force is their per-term forces plus the Pulay sum -- Forces() itself refuses while the XC term has no face,
    // so the test sums MakeTermForces and Tr(W ∂S/∂R) the way MakeForces does, skipping only that zero.  (The Dunlap term's own fixed-D
    // gate is TheDunlapHartreesForceAtFixedDIsItsEnergysFiniteDifference, since the fit's target is the fixed N.)  Each displaced point is a
    // full SCF (MinΔρ/N 1e-9, energy error second order); h = 1e-4 puts the FD's h²/6 |E‴| at ≲ 1e-8, so 2e-6 (relative above 1).
    // _Long: 18 displaced SCFs, ~59 s CPU serial (2026-11-08).
    const Molecule h2o=GetMolecule("H2O");
    const CalcOptions opts{.basis="dzvp", .model=Model::Xalpha, .xalpha=0.0};
    MoleculeCalculation c(MoleculeSetup::Resolve(h2o, opts), Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const auto Ft=ForceFacadeTests::TermForcesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    const auto Et=ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    ASSERT_EQ(Ft.size(), Et.size());
    rvec3vec_t F=ForceFacadeTests::Pulay(c, ForceFacadeTests::EnergyWeighted(c));
    size_t missing=0; bool dunlap=false;
    for (size_t t=0;t<Ft.size();t++)
    {
        if (!Ft[t].missing.empty())
        {
            ++missing;
            EXPECT_EQ(Et[t].E, 0.0) << "only the α = 0 XC term may lack a force: " << Ft[t].missing;
            continue;
        }
        if (Et[t].entries=="Eee") dunlap=true;
        for (size_t A=0;A<F.size();A++) F[A]+=Ft[t].F[A];
    }
    EXPECT_EQ(missing, 1u) << "the XC term";
    EXPECT_TRUE(dunlap) << "the Dunlap Hartree has its force";
    const double h=1e-4;
    rvec3_t sum(0,0,0);
    double largest=0.0;
    for (size_t A=0;A<F.size();A++)
    {
        sum+=F[A];
        for (int x=0;x<3;x++)
        {
            auto E=[&](double s){ MoleculeCalculation m(MoleculeSetup::Resolve(dynamic_cast<const Molecule&>(*MovedH2O(h2o,A,x,s)), opts), Tight());
                                  EXPECT_TRUE(m.IsConverged()); return m.Energy(); };
            const double fd=(E(h)-E(-h))/(2*h);
            EXPECT_NEAR(C(F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd))) << "atom " << A << " axis " << x;
            largest=std::max(largest, std::fabs(fd));
        }
    }
    EXPECT_GT(largest, 1e-2) << "the force is not trivially small (a sign would show)";
    EXPECT_LT(norm(sum), 1e-8) << "translation invariance";
}

TEST(HamiltonianForces, UnderAPseudopotentialEachTermsForceAtFixedDIsItsEnergysFiniteDifference)
{
    // GTH-LDA water, DZVP, closed-form PP (GeometryOptimizationPlan §13.2): the local PP, the KB projectors, the kinetic term, the ion-ion
    // term (valence charges) and the Dunlap Hartree each against the central difference of its own energy at R ± h with the converged D
    // frozen (the Dunlap fit re-solved at each point under the FIXED N, so no λ piece: 2026-11-09).  The XC term has no force yet, and is
    // the only term missing one.  Each is smooth in R: at h = 1e-4 the h²/6 |E‴| is ≲ 1e-8, so 2e-6.
    const Molecule h2o=GetMolecule("H2O");
    MoleculeCalculation c(MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"}), Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const auto Ft=ForceFacadeTests::TermForcesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    const auto Et=ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    ASSERT_EQ(Ft.size(), Et.size());
    std::vector<size_t> gated; size_t missing=0; bool local=false, kb=false;
    for (size_t t=0;t<Ft.size();t++)
    {
        if (!Ft[t].missing.empty()) { ++missing; std::cout<<"[forces] PP water still missing: "<<Ft[t].missing<<std::endl; continue; }
        local |= Et[t].label.rfind("Local pseudopotential",0)==0;
        kb    |= Et[t].label.rfind("Separable (KB)",0)==0;
        gated.push_back(t);
    }
    EXPECT_EQ(missing, 1u) << "only the XC term";
    EXPECT_TRUE(local && kb) << "both PP terms have their force";
    const double h=1e-4;
    for (size_t t : gated)
    {
        rvec3_t sum(0,0,0);
        for (size_t A=0;A<Ft[t].F.size();A++)
        {
            sum+=Ft[t].F[A];
            for (int x=0;x<3;x++)
            {
                auto E=[&](double s){ return ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *MovedH2O(h2o,A,x,s), D)[t].E; };
                const double fd=(E(h)-E(-h))/(2*h);
                EXPECT_NEAR(C(Ft[t].F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd))) << Et[t].label << ", atom " << A << " axis " << x;
            }
        }
        EXPECT_LT(norm(sum), 1e-8) << Et[t].label << ": translation invariance";
    }
}

TEST(HamiltonianForces, TheDunlapHartreesForceAtFixedDIsItsEnergysFiniteDifference)
{
    // All-electron LDA water, DZVP, the default (Dunlap) Hartree: its force at the converged D frozen against the central difference of
    // its own energy entry, the fit re-solved at R ± h.  The fit is charge-constrained to the FIXED electron count N (user ruling
    // 2026-11-08, landed with lane S's 3b/3c), so stationarity in c leaves λ dN/dR = 0 and the fixed-D derivative IS the term's force --
    // the gate this term could not have while the target was Tr(DS), which moved with the basis (2e-4 relative on water then).  Smooth
    // in R: at h = 1e-4 the h²/6 |E‴| is ≲ 1e-8, so 2e-6 (relative above 1).
    const Molecule h2o=GetMolecule("H2O");
    MoleculeCalculation c(MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::LDA}), Tight());
    ASSERT_TRUE(c.IsConverged());
    const auto D=ForceFacadeTests::Carry(c);
    const auto Ft=ForceFacadeTests::TermForcesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    const auto Et=ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
    ASSERT_EQ(Ft.size(), Et.size());
    size_t t=Ft.size();
    for (size_t i=0;i<Ft.size();i++) if (Et[i].entries=="Eee") t=i;
    ASSERT_LT(t, Ft.size()) << "no Hartree (Eee) term";
    ASSERT_TRUE(Ft[t].missing.empty()) << Ft[t].missing;
    const double h=1e-4;
    rvec3_t sum(0,0,0);
    double largest=0.0;
    for (size_t A=0;A<Ft[t].F.size();A++)
    {
        sum+=Ft[t].F[A];
        for (int x=0;x<3;x++)
        {
            auto E=[&](double s){ return ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *MovedH2O(h2o,A,x,s), D)[t].E; };
            const double fd=(E(h)-E(-h))/(2*h);
            EXPECT_NEAR(C(Ft[t].F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd))) << Et[t].label << ", atom " << A << " axis " << x;
            largest=std::max(largest, std::fabs(fd));
        }
    }
    EXPECT_GT(largest, 1e-2) << "the force is not trivially small (a sign would show)";
    EXPECT_LT(norm(sum), 1e-8) << "translation invariance";
}

// ============================================================================================================================
// THE HF FORCE (§13.2; lane B 2026-11-08): the Coulomb and exchange terms through the density's sweep pullbacks, so every HF term has one.
// ============================================================================================================================
TEST(HamiltonianForces, TheHFCoulombAndExchangeForcesAtFixedDAreTheirEnergysFiniteDifference)
{
    // HF water, DZVP, RHF and UHF: the Coulomb (Eee = ½ Tr DJ) and exchange (Exc = ½ Σ_σ Scale Tr D_σ K_σ) forces at the converged D
    // frozen, against the central difference of each one's energy entry at R ± h.  Both are quadratic in D and their matrices are their
    // energies' D-gradients, so nothing is left out.  Smooth in R: at h = 1e-4 the h²/6 |E‴| is ≲ 1e-8, so 2e-6 (relative above 1).
    const Molecule h2o=GetMolecule("H2O");
    for (SpinGroup spin : {SpinGroup::UnPolarized, SpinGroup::Polarized})
    {
        MoleculeCalculation c(MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::HF, .spin=spin}), Tight());
        ASSERT_TRUE(c.IsConverged());
        EXPECT_TRUE(ForceFacadeTests::MissingGradients(c).empty()) << "every HF term has its force";
        const auto D=ForceFacadeTests::Carry(c);
        const auto Ft=ForceFacadeTests::TermForcesAt(c.SizedSetup(), ForceFacadeTests::Structure_(c), D);
        for (const auto& [label, entry] : {std::pair<std::string,std::string>{"Coulomb potential", "Eee"}, {"Hartree-Fock exchange", "Exc"}})
        {
            const Hamiltonian::TermForce* t=nullptr;
            for (const auto& f : Ft) if (f.label.rfind(label,0)==0) t=&f;
            ASSERT_NE(t, nullptr) << "no term labelled " << label;
            ASSERT_TRUE(t->missing.empty()) << t->missing;
            const double h=1e-4;
            rvec3_t sum(0,0,0);
            for (size_t A=0;A<t->F.size();A++)
            {
                sum+=t->F[A];
                for (int x=0;x<3;x++)
                {
                    auto E=[&](double s){ return ForceFacadeTests::TotalEnergyAt(c.SizedSetup(), *MovedH2O(h2o,A,x,s), D)[entry]; };
                    const double fd=(E(h)-E(-h))/(2*h);
                    EXPECT_NEAR(C(t->F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd)))
                        << (spin==SpinGroup::Polarized ? "UHF " : "RHF ") << label << ", atom " << A << " axis " << x;
                }
            }
            EXPECT_LT(norm(sum), 1e-8) << label << ": translation invariance";
        }
    }
}

TEST(HamiltonianForces, TheHFForceIsTheFiniteDifferenceOfTheReconvergedEnergy)
{
    // RHF water, DZVP: every term has its force, so Forces() is the whole thing (the Pulay sum included), against the central difference
    // of the re-converged energy.  Each displaced point is a full SCF (MinΔρ/N 1e-9, energy error second order); h = 1e-4 puts the FD's
    // h²/6 |E‴| at ≲ 1e-8, so 2e-6 (relative above 1).  The stored water sits within ~1e-4 Ha/bohr of its HF/DZVP minimum, where that
    // tolerance could not tell a sign, so one H is pulled 0.15 bohr off it: forces of ~1e-1, checked non-trivial below.
    const Molecule h2o=dynamic_cast<const Molecule&>(*MovedH2O(GetMolecule("H2O"), 1, 0, 0.15));
    const CalcOptions opts{.basis="dzvp", .model=Model::HF};
    MoleculeCalculation c(MoleculeSetup::Resolve(h2o, opts), Tight());
    ASSERT_TRUE(c.IsConverged());
    const rvec3vec_t F=ForceFacadeTests::Forces(c, ForceFacadeTests::Carry(c), ForceFacadeTests::EnergyWeighted(c));
    const double h=1e-4;
    rvec3_t sum(0,0,0);
    double largest=0.0;
    for (size_t A=0;A<F.size();A++)
    {
        sum+=F[A];
        for (int x=0;x<3;x++)
        {
            auto E=[&](double s){ MoleculeCalculation m(MoleculeSetup::Resolve(dynamic_cast<const Molecule&>(*MovedH2O(h2o,A,x,s)), opts), Tight());
                                  EXPECT_TRUE(m.IsConverged()); return m.Energy(); };
            const double fd=(E(h)-E(-h))/(2*h);
            EXPECT_NEAR(C(F[A],x), -fd, 2e-6*std::max(1.0,std::fabs(fd))) << "atom " << A << " axis " << x;
            largest=std::max(largest, std::fabs(fd));
        }
    }
    EXPECT_GT(largest, 1e-2) << "the force is not trivially small (a sign would show)";
    EXPECT_LT(norm(sum), 1e-8) << "translation invariance";
}

// ============================================================================================================================
// A SYMMETRY-BLOCKED RUN'S FORCE (§13.2; lane B 2026-11-09): each SALC block's pullback is the raw AO basis's on its density folded up
// (O D Oᵀ), and the HF sweep folds every block first and pulls back once.
// ============================================================================================================================
namespace
{
//! \a h2o scaled by \a f about its nuclear centroid: a similarity about a point every symmetry element holds, so the point group is kept.
Molecule StretchedH2O(const Molecule& h2o, double f)
{
    std::vector<rvec3_t> R; h2o.ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
    rvec3_t c(0,0,0); for (const auto& r : R) c+=r; c/=double(R.size());
    for (auto& r : R) r=c+f*(r-c);
    std::shared_ptr<Structure> st=h2o.Clone(); st->MoveSites(R);
    return dynamic_cast<const Molecule&>(*st);
}
}

TEST(HamiltonianForces, ASymmetryBlockedRunsForceIsTheUnblockedRunsForce)
{
    // SALC blocking (.symmetry=true) changes how D is stored, not the physics, so at the same converged density the blocked run's force is
    // the unblocked run's -- whose own gates above check it against the finite difference.  Water stretched 8% about its centroid (C2v kept,
    // forces ~1e-1, checked non-trivial).  Every term's force at fixed D is compared, and the whole force (the Pulay sum included) where the
    // model has every term's.  Both SCFs reach MinΔρ/N 1e-9 and a force is first order in the density's error, so 1e-6.  The cases cover
    // the 1E faces (E1), the HF whole-system sweep (RHF, UHF), and under LDA the Dunlap 3-centre and the PP faces; LDA's XC term has no
    // force in either run.
    const Molecule h2o=StretchedH2O(GetMolecule("H2O"), 1.08);
    const std::vector<CalcOptions> cases{{.basis="dzvp", .model=Model::E1}, {.basis="dzvp", .model=Model::HF},
                                         {.basis="dzvp", .model=Model::HF, .spin=SpinGroup::Polarized}, {.basis="dzvp", .model=Model::LDA},
                                         {.basis="dzvp", .model=Model::LDA, .pp="gth-lda-valence"}};
    for (CalcOptions free : cases)
    {
        CalcOptions blocked=free; blocked.symmetry=true;
        const std::string what=std::string(free.model==Model::E1 ? "E1" : free.model==Model::HF ? "HF" : "LDA")
                              +(free.spin==SpinGroup::Polarized ? " polarized" : "")+(free.pp.empty() ? "" : " PP");
        MoleculeCalculation f(MoleculeSetup::Resolve(h2o, free), Tight()), b(MoleculeSetup::Resolve(h2o, blocked), Tight());
        ASSERT_TRUE(f.IsConverged() && b.IsConverged()) << what;
        ASSERT_EQ(ForceFacadeTests::NumBlocks(f), 1u) << what;
        ASSERT_GT(ForceFacadeTests::NumBlocks(b), 1u) << what << ": the run is blocked";
        EXPECT_NEAR(b.Energy(), f.Energy(), 1e-8) << what << ": the same state";
        EXPECT_EQ(ForceFacadeTests::MissingGradients(b), ForceFacadeTests::MissingGradients(f)) << what;

        const auto Ff=ForceFacadeTests::TermForcesAt(f.SizedSetup(), ForceFacadeTests::Structure_(f), ForceFacadeTests::Carry(f));
        const auto Fb=ForceFacadeTests::TermForcesAt(b.SizedSetup(), ForceFacadeTests::Structure_(b), ForceFacadeTests::Carry(b));
        ASSERT_EQ(Fb.size(), Ff.size()) << what;
        double largest=0.0; size_t compared=0;
        for (size_t t=0;t<Ff.size();t++)
        {
            ASSERT_EQ(Fb[t].label, Ff[t].label) << what;
            EXPECT_EQ(Fb[t].missing.empty(), Ff[t].missing.empty()) << what << ", " << Ff[t].label << ": " << Fb[t].missing;
            if (!Fb[t].missing.empty() || !Ff[t].missing.empty()) continue;
            ++compared;
            for (size_t A=0;A<Ff[t].F.size();A++) for (int x=0;x<3;x++)
            {
                EXPECT_NEAR(C(Fb[t].F[A],x), C(Ff[t].F[A],x), 1e-6) << what << ", " << Ff[t].label << ", atom " << A << " axis " << x;
                largest=std::max(largest, std::fabs(C(Ff[t].F[A],x)));
            }
        }
        EXPECT_GE(compared, 3u) << what;
        EXPECT_GT(largest, 1e-2) << what << ": the terms' forces are not trivially small";
        if (!ForceFacadeTests::MissingGradients(f).empty()) continue;
        const rvec3vec_t Wf=ForceFacadeTests::Forces(f, ForceFacadeTests::Carry(f), ForceFacadeTests::EnergyWeighted(f));
        const rvec3vec_t Wb=ForceFacadeTests::Forces(b, ForceFacadeTests::Carry(b), ForceFacadeTests::EnergyWeighted(b));
        double total=0.0;
        for (size_t A=0;A<Wf.size();A++) for (int x=0;x<3;x++)
        {
            EXPECT_NEAR(C(Wb[A],x), C(Wf[A],x), 1e-6) << what << ", the whole force, atom " << A << " axis " << x;
            total=std::max(total, std::fabs(C(Wf[A],x)));
        }
        EXPECT_GT(total, 1e-2) << what << ": the force is not trivially small";
    }
}

// ============================================================================================================================
// THE HAND-GRADIENT AUDIT (lane B 2026-11-09): the converged state's fields' Gradient(r) -- ρ (the composite of irrep densities, each
// Σ D ∇φ φ), m = ρ↑−ρ↓ (the spin density), and every MO (Σ c ∇φ) -- against the central difference of their own values, RHF and UHF
// water, free and SALC-blocked.  Step h = 1e-5 and tolerance as M_GeometryDerivatives GATE 2 (1e-7, relative above 1); the largest
// derivative probed must be 0.1 or more (the spin density's 0.01: water's UHF m is small), so a sign would show.
// ============================================================================================================================
TEST(HandGradient, TheConvergedStatesFieldsGradientsAreTheFiniteDifferenceOfTheirValues)
{
    const Molecule h2o=GetMolecule("H2O");
    std::vector<rvec3_t> pts; h2o.ForEachSite([&](int, const rvec3_t& r, bool){ pts.push_back(r+rvec3_t(0.23,-0.17,0.31)); });
    pts.push_back(rvec3_t(1.4,-0.9,2.2));
    auto check=[&](const ScalarFunction<double>& f, const std::string& who, double floor)
    {
        const double h=1e-5;
        double largest=0.0;
        for (const rvec3_t& r : pts)
        {
            const rvec3_t g=f.Gradient(r);
            for (int k=0;k<3;k++)
            {
                rvec3_t p=r, m=r; (k==0 ? p.x : k==1 ? p.y : p.z)+=h; (k==0 ? m.x : k==1 ? m.y : m.z)-=h;
                const double fd=(f(p)-f(m))/(2*h);
                EXPECT_NEAR(C(g,k), fd, 1e-7*std::max(1.0,std::fabs(fd))) << who << ": axis " << k << " at " << r;
                largest=std::max(largest, std::fabs(fd));
            }
        }
        EXPECT_GT(largest, floor) << who << ": the derivatives probed are not trivially small";
    };
    for (SpinGroup spin : {SpinGroup::UnPolarized, SpinGroup::Polarized})
        for (bool sym : {false, true})
        {
            const std::string what=std::string(spin==SpinGroup::Polarized ? "UHF" : "RHF")+(sym ? " SALC" : " free");
            MoleculeCalculation c(MoleculeSetup::Resolve(h2o, {.basis="dzvp", .model=Model::HF, .spin=spin, .multiplicity=spin==SpinGroup::Polarized ? 3 : 0, .symmetry=sym}), Tight());
            ASSERT_TRUE(c.IsConverged()) << what;
            check(c.Density(), what+" density", 0.1);
            if (spin==SpinGroup::Polarized) check(c.GetSpinDensity(), what+" spin density", 0.01);
            size_t n=0;
            for (const auto& L : c.GetLevels()) if (n++<8) check(*L.phi, what+" MO "+L.label, 0.1);
        }
}

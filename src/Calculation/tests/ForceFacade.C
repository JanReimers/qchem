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
// it sums to zero over the atoms, and its AD gradient in D IS the term's H (pin 32).  Built without Enzyme, the gradients refuse.
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
    if (!AD::Available()) { ASSERT_FALSE(g) << "a build without Enzyme refuses"; std::cout<<"[AD spike] "<<g.Error().why<<std::endl; return; }
    ASSERT_TRUE(g);
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
    ASSERT_TRUE(g2 && gD);
    for (size_t i=0;i<R.size();i++) EXPECT_EQ((*g)[i], (*g2)[i]) << "the gradient is deterministic (block sums in block order)";
    // THE PER-TERM FINITE DIFFERENCE: the library's delta XC term at R +- h with D frozen (basis centres, mesh points and partition move).
    const double h=1e-4;
    rvec3_t sum(0,0,0);
    for (size_t A=0;A<d.natom;++A)
    {
        sum+=rvec3_t((*g)[3*A], (*g)[3*A+1], (*g)[3*A+2]);
        for (int x=0;x<3;++x)
        {
            auto at=[&](double s){ std::vector<rvec3_t> Rv; for (size_t B=0;B<d.natom;B++) Rv.push_back(rvec3_t(R[3*B],R[3*B+1],R[3*B+2]));
                                    (x==0 ? Rv[A].x : x==1 ? Rv[A].y : Rv[A].z)+=s; std::shared_ptr<Structure> st=h2o.Clone(); st->MoveSites(Rv);
                                    return Entry(ForceFacadeTests::TermEnergiesAt(c.SizedSetup(), *st, D), "Exc"); };
            const double fd=(at(h)-at(-h))/(2*h);
            std::cout<<"[AD spike] dExc/dR atom "<<A<<" axis "<<x<<": AD "<<std::setprecision(10)<<(*g)[3*A+x]<<"  FD "<<fd<<"  diff "<<(*g)[3*A+x]-fd<<std::endl;
            EXPECT_NEAR((*g)[3*A+x], fd, 2e-6*std::max(1.0, std::fabs(fd))) << "dExc/dR atom " << A << " axis " << x;
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
    if (!AD::Available()) { ASSERT_FALSE(g) << "a build without Enzyme refuses"; return; }
    ASSERT_TRUE(g);
    const hmat_t<double> H=ForceFacadeTests::TermMatrix(c, D, "Exc");
    ASSERT_EQ(H.rows(), d.nbf);
    double worst=0.0, largest=0.0;
    for (size_t i=0;i<d.nbf;i++) for (size_t j=0;j<d.nbf;j++)
    {
        worst=std::max(worst, std::fabs((*g)(i,j)-H(i,j))); largest=std::max(largest, std::fabs(H(i,j)));
        EXPECT_NEAR((*g)(i,j), H(i,j), 1e-11) << "dE/dD_" << i << j << " against the term's H";
    }
    std::cout<<"[AD spike] pin 32: max|dE/dD - H| = "<<worst<<" (max|H| "<<largest<<")"<<std::endl;
}

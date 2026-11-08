// File: Calculation/tests/ForceFacade.C  The facade-side pieces of the analytic force (doc/GeometryOptimizationPlan.md §4-§6, §13),
// with no Hamiltonian-hub change: the energy-weighted density W (the Pulay term's input), and the PER-TERM FINITE DIFFERENCE AT FIXED D
// -- one term's energy at R ± h with the AO density frozen, no SCF -- gated on the nuclear attraction, whose analytic gradient libcint
// gives (int1e_ipnuc + int1e_iprinv, PG_LibCint::NuclearGradient), and on the kinetic term (int1e_ipkin).
#include <gtest/gtest.h>
#include <cmath>
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

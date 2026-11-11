// File: Calculation/tests/MolecularRelaxation.C  REAL MOLECULAR RELAXATIONS (doc/OpenWork.md §2 "Real relaxations end to end"): the
// geometry optimizer over a MoleculeSession, SCF-driven, analytic forces.  HF/DZVP, because every HF term has its analytic force today
// (LDA's XC term does not: MissingGradients names it).  The oracle is PySCF at the SAME model: the Cartesian DGauss DZVP basis read from
// BasisSetData/dzvp.bsd, RHF, analytic gradients, scipy BFGS to max|g| < 1e-7 Ha/bohr (~/Code/pyscf-runs/Relax/hf_dzvp_relax.py, its
// output beside it).  The starting geometries are experimental (H2O) and a torsion 54 degrees off the minimum (H2O2).  NH3 (the same
// recipe, 4 steps, E and geometry on PySCF's to 1e-10 Ha / 3e-6 bohr) is banked in doc/GeometryOptimizationPlan.md §15, not gated.
#include <gtest/gtest.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
import qchem.MoleculeCalculation;
import qchem.GeometryOptimizer;
import qchem.StructureData;
import qchem.Structure;
import qchem.SCFParams;
import qchem.Types;
import qchem.Reporting;   // report::FormatGuard

using namespace qchem;

namespace
{
SCFParams Tight() {SCFParams p; p.NMaxIter=100; p.MinΔρPerElectron=1e-9; return p;}

struct Atom_ {int Z; rvec3_t R;};
Molecule Make(const std::vector<Atom_>& atoms)
{
    Molecule m;
    for (const auto& a : atoms) m.Insert(Site{.Z=a.Z, .R=a.R});
    return m;
}

double Distance(const std::vector<rvec3_t>& R, size_t i, size_t j) {return norm(R[i]-R[j]);}
double Angle(const std::vector<rvec3_t>& R, size_t i, size_t j, size_t k)   //!< degrees, at j
{
    const rvec3_t a=R[i]-R[j], b=R[k]-R[j];
    return std::acos((a*b)/(norm(a)*norm(b)))*180.0/M_PI;
}
double Dihedral(const std::vector<rvec3_t>& R, size_t i, size_t j, size_t k, size_t l)   //!< degrees, i-j-k-l
{
    rvec3_t b1=R[k]-R[j]; b1/=norm(b1);
    const rvec3_t b0=R[i]-R[j], b2=R[l]-R[k];
    const rvec3_t v=b0-(b0*b1)*b1, w=b2-(b2*b1)*b1;
    return std::atan2(Cross(b1, v)*w, v*w)*180.0/M_PI;
}

//! Relax \a m from its given positions; print the history; hand back the result (asserting it converged).
OptimizedGeometry Relax(const Molecule& m, OptimizerOptions::Algorithm algorithm, const std::string& what)
{
    MoleculeSession pes(MoleculeSetup::Resolve(m, {.basis="dzvp"}), Tight());
    CartesianMap map(pes.Reference());
    OptimizerOptions o=OptimizerOptions::FromEps(1e-6);
    o.algorithm=algorithm;
    GeometryOptimizer opt(pes, map, o);
    auto r=opt.Run();
    EXPECT_TRUE(bool(r)) << what << ": " << (r ? std::string() : r.Error().details);
    if (!r) return {};
    report::FormatGuard keep(std::cout);
    std::cout<<"[relax] "<<what<<(algorithm==OptimizerOptions::Algorithm::L_BFGS ? " L-BFGS" : " BFGS-RFO")<<": "<<r->Steps()
             <<" steps, "<<r->history.size()<<" evaluations"<<std::endl;
    for (const auto& h : r->history)
        std::cout<<"[relax]   "<<std::setw(3)<<h.index<<"  E "<<std::setprecision(10)<<std::fixed<<h.value<<std::defaultfloat
                 <<"  max|g| "<<std::setprecision(3)<<std::scientific<<h.maxGradient<<"  |s| "<<h.stepLength
                 <<std::defaultfloat<<(h.accepted ? "" : "  rejected")<<std::endl;
    return r.TakeValue();
}
} // namespace

// H2O, HF/DZVP, from the experimental geometry.  PySCF: E -76.0232244736, r(OH) 1.788807 bohr, HOH 106.4775 degrees.
TEST(MolecularRelaxation, HFWaterReachesPySCFsMinimum)
{
    const Molecule h2o=StructureData::GetMolecule("H2O");
    for (auto alg : {OptimizerOptions::Algorithm::BFGS_RFO, OptimizerOptions::Algorithm::L_BFGS})
    {
        const auto r=Relax(h2o, alg, "H2O");
        const auto& R=r.geometry.R;
        std::cout<<"[relax]   E "<<std::setprecision(12)<<r.value<<" (PySCF -76.0232244736)  r(OH) "<<Distance(R,0,1)<<" "<<Distance(R,0,2)
                 <<" (1.788807)  HOH "<<Angle(R,1,0,2)<<" (106.4775)"<<std::endl;
        EXPECT_NEAR(r.value, -76.0232244736, 1e-7);
        EXPECT_NEAR(Distance(R,0,1), 1.788807, 1e-4);
        EXPECT_NEAR(Distance(R,0,2), 1.788807, 1e-4);
        EXPECT_NEAR(Angle(R,1,0,2), 106.4775, 0.05);
    }
}

// H2O2, HF/DZVP, the SOFT TORSION: started at HOOH 60 degrees (O-O 1.40 A, OOH 102), relaxed to PySCF's E -150.7819917055, r(OO) 2.636291,
// r(OH) 1.792327 bohr, OOH 102.7434, HOOH 113.7064 degrees.  The molecular family's ONE end-to-end relaxation (user 2026-11-10, the test
// budget), 13 points.  It was `_Long` at ~150 s, ~80 % of it the MnD ERI4 rebuilt per point; since the unique-quartet ERI4 (2026-11-11,
// OpenWork §4) it is ~32 s, so it runs in the short tier.
TEST(MolecularRelaxation, HFPeroxideUnwindsItsTorsionToPySCFsMinimum)
{
    const Molecule h2o2=Make({{8,{-1.322808,0,0}}, {8,{1.322808,0,0}}, {1,{-0.949557,1.520749,0.878005}}, {1,{0.949557,-1.520749,0.878005}}});
    const auto r=Relax(h2o2, OptimizerOptions::Algorithm::BFGS_RFO, "H2O2");
    const auto& R=r.geometry.R;
    std::cout<<"[relax]   E "<<std::setprecision(12)<<r.value<<" (PySCF -150.7819917055)  r(OO) "<<Distance(R,0,1)<<" (2.636291)  r(OH) "<<Distance(R,0,2)
             <<" (1.792327)  OOH "<<Angle(R,1,0,2)<<" (102.7434)  HOOH "<<Dihedral(R,2,0,1,3)<<" (113.7064)"<<std::endl;
    EXPECT_NEAR(r.value, -150.7819917055, 1e-7);
    EXPECT_NEAR(Distance(R,0,1), 2.636291, 1e-4);
    EXPECT_NEAR(Distance(R,0,2), 1.792327, 1e-4);
    EXPECT_NEAR(std::fabs(Dihedral(R,2,0,1,3)), 113.7064, 0.5);
}

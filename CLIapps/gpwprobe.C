// File: CLIapps/gpwprobe.C  THE GPW CAMPAIGN INSTRUMENTS -- the hand-run ladders, sweeps and probes that used to
// sit in IntegrationTests/GPW_SCF_UT.C as DISABLED_ tests (doc/TestSuitePlan.md §8, verdict P, 2026-09-15).
//
// An instrument is something a human runs with knobs and READS; a test is something ctest runs and JUDGES.
// The two had grown together in one file, which is why nobody could say which of 62 GPW_SCF entries were
// gates.  Every sub-command here drives qchem::SolidCalculation through the same harness the gates use
// (IntegrationTests/GPW/Harness.C), so a probe and the gate it informs build the same cell the same way.
// The env-var knobs are kept VERBATIM: doc/Benchmark.md's recipes are written with them, and a recipe you
// have to reconstruct is a recipe you get wrong (the copy-the-command rule).
//
//   gpwprobe ladder                SI_LADDER=n1,n2,n3  SI_XC=becke      the Si supercell scaling ladder (Γ)
//   gpwprobe ksweep                GPW_KSHIFT=s                         E(k) at ONE k = s(1,1,1) on Si
//   gpwprobe naf-smear             NAFGDM_*                             NaF IBZ, DIIS-smear then cold GDM
//   gpwprobe becke-ladder SYSTEM   SYSTEM = si | naf | mn | al          the V2.6 Becke (nR, degree) ladder
//   (mno / nio: RETIRED 2026-10-04 -> `rundeck decks/MnO_AFM2_free.json`, `decks/NiO_AFM2_acbn0.json`)
//
// Each sub-command prints its findings and the CHECKS the retired test asserted (PASS/FAIL lines); the exit
// code is the number of failed checks.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <memory>
#include <cmath>
#include <complex>
#include <stdexcept>

import qchem.Tests.GPW_Harness;      // the recipe vocabulary + the instruments shared with the gates
import qchem.UnitCell;               // Supercell, BravaisCell (the MnO cell's named construction)
import qchem.Matrix3D;
import qchem.ChargeDensity;          // cSpinResolved_CD
import qchem.ChargeDensity.FourierDensity;
import qchem.ScalarFunction;
import qchem.Symmetry.Spin;
import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type
import qchem.Blaze;                    // the complex arithmetic on the ΔG_Map entries
import qchem.Types;

using namespace qchem;
using namespace qchem::tests::gpw;

namespace
{
int nFailed=0;
void Check(bool ok, const std::string& what)
{
    std::cout << (ok ? "[PASS] " : "[FAIL] ") << what << std::endl;
    if (!ok) ++nFailed;
}
void CheckNear(double x, double ref, double tol, const std::string& what)
{
    std::ostringstream s; s << what << ": " << std::setprecision(10) << x << " vs " << ref << " (tol " << tol << ")";
    Check(std::fabs(x-ref)<=tol, s.str());
}
double Envd(const char* n, double d) { const char* s=std::getenv(n); return s ? std::atof(s) : d; }
int    Envi(const char* n, int    d) { const char* s=std::getenv(n); return s ? std::atoi(s) : d; }

//========================================================================================================
// ★ THE SUPERCELL SCALING LADDER (doc/ParallelAndOraclePlan.md 2.1).  Every threading number we own was
// measured on a 4-atom cell, and a small cell starves threads on both sides of the CP2K comparison -- so
// the question this answers is whether our speedup is a property of the CODE or of the development cell.
// ONE RUNG PER INVOCATION (`SI_LADDER=n1,n2,n3`, default 1,1,1): peak RSS is a PROCESS watermark, so a
// single process walking the whole ladder would report one number for the largest rung.
// ★★ AND IT IS A CORRECTNESS GATE FOR FREE, WHICH IS WHY THE LADDER IS AT Γ.  A Γ-only calculation on an
// N1xN2xN3 SUPERCELL is band-folding-equivalent to an N1xN2xN3 k-MESH on the primitive cell, so each rung
// must reproduce the k-mesh total the suite banks PER PRIMITIVE CELL:
//   1x1x1 -> -7.11506  (GPW_Si.Γ_Imp_CP2K)   2x1x1 -> -7.45294  (GPW_Si.k211_Imp_Anchor, KP-0 re-pin)
//   2x2x2 -> -7.77846  (GPW_Si.k222_CP2K)     (2x2x1 has no banked counterpart: timing only)
// SI_XC=becke forces the atom-centred mesh -- the path that exercises SiteStabilizer in the supercell
// setting (the §6a W2b site-adapted angular sets).  Setting cellKind ALONE is a trap: ask for the RECIPE
// (BeckeXCParams), which also makes QCHEM_BECKE_L / QCHEM_BECKE_NR live.  The banked anchors are UNIFORM-mesh
// numbers; a Becke rung is a different quadrature (~75 mHa away at the coarse default) and is not compared.
//========================================================================================================
int Ladder()
{
    const Material si=Materials::Get("Si_diamond");
    ivec3_t n(1,1,1);
    if (const char* e=std::getenv("SI_LADDER")) { int x=1,y=1,z=1; sscanf(e,"%d,%d,%d",&x,&y,&z); n=ivec3_t(x,y,z); }
    UnitCell cell=Supercell(*si.cell, n);
    const size_t nAtom=cell.GetNumAtoms();
    const int    Nelec=4*int(nAtom);                 // Si Zion=4 via the PP
    const size_t nPrim=nAtom/2;

    const bool becke = std::getenv("SI_XC") && std::string(std::getenv("SI_XC"))=="becke";
    std::ostringstream label; label<<"Si supercell "<<n.x<<"x"<<n.y<<"x"<<n.z<<" ("<<nAtom<<" atoms) Gamma";
    const Lattice_3D lat(cell, ivec3_t(1,1,1));      // Γ ONLY -- the folding equivalence above

    SolidCalcOptions o;
    o.label=label.str(); o.Nelec=Nelec; o.species=si.species;
    o.densityEcut=20.0; o.imposeSymmetry=true;
    o.seed=ChargeDensity::SeedStrategy::Uniform;
    if (becke) { o.xcMesh=qcMesh::BeckeXCParams(-1,-1.0,-1); o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke; }
    SCFParams par=ProductionGates();
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    SolidCalculation calc(lat, MakeBasisSR(cell), o, par);
    auto R=calc.Result();
    const double E=calc.LastIterateTerms().GetTotalEnergy(), ePerPrim=E/double(nPrim);
    std::cout<<"[ladder] "<<n.x<<"x"<<n.y<<"x"<<n.z<<"  atoms="<<nAtom<<"  Etot="<<std::setprecision(10)<<E
             <<"  E/primitive="<<ePerPrim<<"  iters="<<calc.IterationCount()<<std::endl;
    Check(bool(R), "converged" + (R ? std::string() : ": "+R.Error().details));
    CheckNear(calc.LastIterateCharge(), double(Nelec), 1e-5, "charge");
    if (becke) return nFailed;
    if (n==ivec3_t(1,1,1)) CheckNear(ePerPrim, -7.11506, 5e-3, "E/primitive vs the Γ anchor");
    if (n==ivec3_t(2,1,1)) CheckNear(ePerPrim, -7.45294, 8e-3, "E/primitive vs the 2x1x1 k-mesh anchor");
    if (n==ivec3_t(2,2,2)) CheckNear(ePerPrim, -7.77846, 1e-2, "E/primitive vs the 2x2x2 k-mesh anchor");
    return nFailed;
}

//========================================================================================================
// SINGLE-K SWEEP: E(k) along the cell diagonal, ONE k-point per run, no weights, no symmetry, no IBZ.  The
// mesh is 1x1x1, so k = kShift EXACTLY: run k = s(1,1,1) for s in [0, 1/2].  s=0 (Γ) and s=1/2 (L) are TRIM
// -- real blocks; everything strictly between is genuinely complex and exercised by nothing else.  E(k) is
// not a physical total (one k is not a BZ average) but it IS smooth in k, so a KINK on entering the complex
// region localizes a defect to a single k-block's operators.
//   for s in 0 0.125 0.25 0.375 0.5; do GPW_KSHIFT=$s gpwprobe ksweep; done
//========================================================================================================
int KSweep()
{
    const Material si=Materials::Get("Si_diamond");
    const Lattice_3D lat=LatticeOf(si);              // ONE k-point, at exactly kShift
    const double s=Envd("GPW_KSHIFT", 0.25);
    const bool trim=(s==0.0 || s==0.5);
    SolidCalcOptions o=OptionsFor(si, "Si single-k s="+std::to_string(s));
    o.densityEcut=20.0; o.imposeSymmetry=true; o.kShift=rvec3_t(s,s,s);
    SCFParams par=TightGates(60);
    EnvOverrides(o, par);
    GpwReport report("Si "+o.label, par.Verbose);
    SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
    const EnergyBreakdown E=calc.LastIterateTerms();
    std::cout << "[k sweep] s=" << s << (trim ? "  TRIM" : "  complex")
              << "  Etot=" << std::setprecision(10) << E.GetTotalEnergy() << std::setprecision(6)
              << "  charge=" << calc.LastIterateCharge()
              << "  (Ekin=" << E["Kinetic"] << " Een=" << E["Een"] << " Eee=" << E["Eee"] << " Exc=" << E["Exc"] << ")" << std::endl;
    CheckNear(calc.LastIterateCharge(), 8.0, 1e-6, "charge (the one thing that must hold at EVERY k)");
    return nFailed;
}

//========================================================================================================
// NaF IBZ, GDM x SMEARING PROBE (doc/SymmetryUpgradePlan.md §6b rounds 1-3): a DIIS-smeared hot stage
// then a COLD GDM stage on the imposed 2x2x2 NaF -- is the imposed E[ρ] variational (GDM strictly
// descends)?  Verdict banked (round 3: yes on SR2; the SR uphill walk is GDM x diffuse near-null modes).
//   NAF_ECUT, NAFGDM_IMPOSE=0/1, NAFGDM_NMAX, NAFGDM_MOM=0/1, NAFGDM_PIVOT=0/1, NAFGDM_L, NAFGDM_KT, NAFGDM_VERBOSE
//========================================================================================================
int NaFSmear()
{
    const Material naf=Materials::Get("NaF_rocksalt");
    const Lattice_3D lat=LatticeOf(naf, ivec3_t(2,2,2));
    SolidCalcOptions o=OptionsFor(naf, "NaF IBZ GDM-smear probe");
    o.densityEcut=Envd("NAF_ECUT",-1.0);                // AUTO = C·αmax = 80 (the anchor config)
    o.imposeSymmetry=Envd("NAFGDM_IMPOSE",1.0)!=0.0;
    o.seed=ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    if (Envd("NAFGDM_PIVOT",1.0)==0.0) { o.ortho=qchem::Cholesky; o.orthoTol=0.0; }   // plain-Cholesky control
    o.xcMesh=qcMesh::BeckeXCParams(20,2,int(Envd("NAFGDM_L",11.0)));
    o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;
    SCFParams base=Gates((size_t)Envd("NAFGDM_NMAX",100.0), 1e-5, 1e-9);   // deep enough to expose a stall vs a clean floor
    base.StartingRelaxRo=0.25; base.KerkerG0=1.0;
    base.UseMOM=Envd("NAFGDM_MOM",1.0)!=0.0; base.MOMStartIter=10;         // NAFGDM_MOM=0: is MOM x GDM the fight?
    base.MergeTol=SCFParams{}.MergeTol;
    base.Verbose=Envd("NAFGDM_VERBOSE",0.0)!=0.0;
    const double kT=Envd("NAFGDM_KT",0.01);
    SCFParams hot=base;  hot.SmearingkT=kT;  hot.StopOnAccelExhausted=true;
    SCFParams cold=base; cold.SmearingkT=0.0;                              // GDM does not smear: its stage runs cold
    const std::vector<SCFStage> schedule={{hot, SCFAccelerators::Type::DIIS}, {cold, SCFAccelerators::Type::GDM}};
    GpwReport report("NaF "+o.label, base.Verbose);
    SolidCalculation calc(lat, MakeBasisLowQ(*naf.cell, BasisSetData::VALENCE_LOWQ_SR), o, schedule);
    auto R=calc.Result();
    std::cout << "[NaF GDM-smear probe] impose="<<o.imposeSymmetry<<" kT(stage1)="<<kT
              << " E="<<std::setprecision(10)<<calc.LastIterateTerms().GetTotalEnergy()
              << " converged="<<bool(R)<<" iters(final)="<<calc.IterationCount()<<std::endl;
    CheckNear(calc.LastIterateCharge(), 8.0, 1e-6, "charge");
    return nFailed;
}

//========================================================================================================
// THE V2.6 BECKE RECIPE LADDER: converge one system, freeze its density, and score a ladder of Becke
// (nR, degree) rules against a strict-refinement reference AND against the uniform route -- the yardstick
// (user 2026-09-06) is the smallest rung whose max|dVxc| matches what the cheap uniform route delivers.
// Four bonding characters: covalent Si, ionic NaF, the open-d Mn sextet box, the nearly-free-electron Al
// metal (the worst case: charge ON the fuzzy-Voronoi partition surface -- it backed the degree-17 flip out).
//========================================================================================================
int BeckeRecipeLadder(const std::string& system)
{
    if (system=="si")
    {
        const Material si=Materials::Get("Si_diamond");
        const Lattice_3D lat=LatticeOf(si);
        SolidCalcOptions o=OptionsFor(si, "Si V2.6 ladder");
        o.densityEcut=60.0;
        o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;   // converge on the uniform route: the probe density
        o.imposeSymmetry=false;                            // free mesh: the ladder measures the RULE, not the fold
        SCFParams par=ProductionGates(); par.MergeTol=SCFParams{}.MergeTol;
        GpwReport report("Si "+o.label, false);
        SolidCalculation calc(lat, MakeBasisSR(*si.cell), o, par);
        Check(bool(calc.Result()), "converged");
        BeckeLadder(Handles(calc), lat.GetStructure(), "Si-covalent");
    }
    else if (system=="naf")
    {
        const Material naf=Materials::Get("NaF_rocksalt");
        const Lattice_3D lat=LatticeOf(naf);
        SolidCalcOptions o=NaFOptions(naf, "NaF V2.6 ladder");
        o.imposeSymmetry=false;
        GpwReport report("NaF "+o.label, false);
        SolidCalculation calc(lat, MakeBasisNaFSR2(*naf.cell), o, NaFGates());
        Check(bool(calc.Result()), "converged");
        BeckeLadder(Handles(calc), lat.GetStructure(), "NaF-ionic");
    }
    else if (system=="mn")
    {
        // NOTE the probe pair is UNPOLARIZED while the run is spin-native: the ladder measures the GRID's
        // ability to integrate a real, strongly aspherical density, not this run's energy.
        const Material box=Materials::Get("Mn_box16");
        const Lattice_3D lat=LatticeOf(box);
        SolidCalcOptions o=MnBoxOptions(box, "Mn V2.6 ladder");
        SCFParams par=Gates(40, 1e-5, 1e30); par.SmearingkT=5e-3;
        std::shared_ptr<const BasisSet::Real_BS> mnbasis(
            BasisSet::Gaussian::Factory(BasisSetData::VALENCE_LOWQ_SR, box.cell.get(), BasisSet::Gaussian::Engine::MnD,
                                        BasisSet::Gaussian::Angular::Cartesian));
        GpwReport report("Mn "+o.label, false);
        SolidCalculation calc(lat, mnbasis, o, par);
        Check(bool(calc.Result()), "converged");
        BeckeLadder(Handles(calc), lat.GetStructure(), "Mn-openshell-d");
    }
    else if (system=="al")
    {
        const Material al=Materials::Get("Al_fcc");
        const Lattice_3D lat=LatticeOf(al);
        SolidCalcOptions o=AlOptions(al, "Al V2.6 ladder");
        o.imposeSymmetry=false;
        SCFParams par=AlGates(); par.SmearingkT=0.02;      // smeared: a converged, non-rotating density to freeze
        GpwReport report("Al "+o.label, false);
        SolidCalculation calc(lat, MakeBasisLowQ(*al.cell, BasisSetData::VALENCE_LOWQ_SR), o, par);
        Check(bool(calc.Result()), "converged");
        BeckeLadder(Handles(calc), lat.GetStructure(), "Al-metal");
    }
    else throw std::runtime_error("becke-ladder: SYSTEM must be si | naf | mn | al");
    return nFailed;
}

//========================================================================================================
// `gpwprobe mno` and `gpwprobe nio` WERE HERE (the rocksalt AFM-II campaign run with ~57 MNO_*/NIO_* environment knobs) -- RETIRED 2026-10-04,
// D-ENV step 6d.3.  A run is now an input DECK:  rundeck decks/MnO_AFM2_free.json  /  decks/NiO_AFM2_acbn0.json  [--set key=value]...  and every
// knob has a deck key (doc/Records/EnvKnobInventory.md §10.2).  Its PASS/FAIL checks live in the integration tests, per the ruling that an
// important gate is a hard-coded test: the AFM-II anchor (GPW_MnO.Γ_Shub_Pol_Smear_Anchor_Long), the ordering gate
// (GPW_MnO.DISABLED_Γ_Shub_Pol_Smear_Ordering_Long -- currently a FALSE claim, doc/OpenWork.md §3), the +U oracle gate (…_U_…_CP2K_Long).
// The probe's m(r) point sample and the |m~(q_AFM)| check were dropped with it (D-POINTPROBE: the integrated site moment replaces them).
//========================================================================================================

//========================================================================================================
// GATE 1 (doc/HubbardUPlan.md track B, B2): does a CHOSEN collinear order on the spinel Mn sublattice --
// the 16d Wyckoff site, i.e. the PYROCHLORE lattice, geometrically frustrated (S3 risk 1) -- survive a
// change in U, consistently and reproducibly?  This is NOT a ground-state search: the NiO lesson (trap 3,
// doc/HubbardUPlan.md S1) is that seeding cannot buy an order the SCF fixed point does not have, so what
// matters is whether the order SURVIVES, not whether it is the true minimum.  Materials come from
// materials.json (LiMn2O4_spinel, MnO2_lambda_spinel), ferromagnetically seeded there already; every Mn
// site gets the SAME Hubbard U (site, l=2).
//   gpwprobe gate1 <material> [U_eV]     GATE1_KT=kT  GATE1_NMAX=n  GATE1_MULT=2S+1 (override the
//                                        d-count-derived default)  GATE1_ECUT
//========================================================================================================
int Gate1(const std::string& materialName, double U_eV)
{
    const Material mat=Materials::Get(materialName);
    std::vector<size_t> mnSites;
    { size_t idx=0; mat.cell->ForEachSite([&](int z, const rvec3_t&, bool){ if (z==25) mnSites.push_back(idx); idx++; }); }
    if (mnSites.empty()) throw std::runtime_error("gate1: no Mn sites in '"+materialName+"'");

    // The formal Mn oxidation state fixes the unambiguous cases (lambda-MnO2: Mn4+, d3, S=3/2) and forces
    // an EXPLICIT approximation on the mixed-valence one (LiMn2O4: formally Mn3.5+ on four CRYSTALLOGRAPHICALLY
    // EQUIVALENT sites -- charge disproportionation is not modelled here).  Ferromagnetic alignment (every
    // Mn spin the same sign, matching materials.json's decoration) is gate 1's own SIMPLE CHOICE, not a
    // ground-state claim.
    const bool isLambda = materialName.find("lambda")!=std::string::npos;
    const int unpairedPerMn = isLambda ? 3 : 4;                     // Mn4+ d3 (S=3/2) : Mn3+ d4 (S=2)
    const int defaultMult = int(mnSites.size())*unpairedPerMn + 1;  // FM: total S = mnSites*unpaired/2
    const int mult = Envi("GATE1_MULT", defaultMult);

    std::ostringstream label; label<<materialName<<" gate1 U="<<U_eV<<"eV mult="<<mult;
    SolidCalcOptions o=OptionsFor(mat, label.str());
    o.multiplicity=mult;
    o.seed=ChargeDensity::SeedStrategy::IonicSAD;
    o.ortho=qchem::CholeskyPivoted; o.orthoTol=1e-4;
    o.cutoffFactor=Envd("GATE1_CUTOFF_FACTOR",2.0);
    o.densityEcut =Envd("GATE1_ECUT",-1.0);
    for (size_t s : mnSites) o.hubbard.push_back(HubbardU(s, 2, U_eV));

    SCFParams par=ProductionGates();
    par.NMaxIter=size_t(Envd("GATE1_NMAX",80.0));
    par.SmearingkT=Envd("GATE1_KT",0.005);
    par.UseMOM=true; par.MOMStartIter=10;
    EnvOverrides(o, par);

    GpwReport report(label.str(), true);
    SolidCalculation calc(LatticeOf(mat), MakeBasisLowQ(*mat.cell, BasisSetData::VALENCE_LOWQ_VA), o, par);
    auto R=calc.Result();
    std::cout << "[" << label.str() << "] " << (R ? "CONVERGED" : "NOT converged")
              << "  Etot=" << std::setprecision(10) << calc.LastIterateTerms().GetTotalEnergy()
              << std::setprecision(6) << "  " << calc.Diagnostics().Summary() << std::endl;
    return nFailed;
}

void Usage()
{
    std::cout << "gpwprobe ladder | ksweep | naf-smear | becke-ladder {si|naf|mn|al} | gate1 <material> [U_eV]\n  (mno / nio are retired: use rundeck with decks/)\n"
                 "  (the env-var knobs are documented at the top of CLIapps/gpwprobe.C and in doc/Benchmark.md)\n";
}
} // anonymous

int main(int argc, char** argv)
{
    if (argc<2) { Usage(); return 2; }
    const std::string cmd=argv[1];
    try
    {
        if (cmd=="ladder")       return Ladder();
        if (cmd=="ksweep")       return KSweep();
        if (cmd=="naf-smear")    return NaFSmear();
        if (cmd=="becke-ladder") { if (argc<3) { Usage(); return 2; } return BeckeRecipeLadder(argv[2]); }
        if (cmd=="mno" || cmd=="nio")
        { std::cerr << "gpwprobe " << cmd << " is RETIRED (D-ENV 6d.3): use  rundeck decks/" << (cmd=="mno" ? "MnO_AFM2_free.json" : "NiO_AFM2_acbn0.json")
                    << "  (every MNO_*/NIO_* knob is a deck key; see doc/Records/EnvKnobInventory.md §10.2)" << std::endl; return 2; }
        if (cmd=="gate1")        { if (argc<3) { Usage(); return 2; }
                                    return Gate1(argv[2], argc>3 ? std::atof(argv[3]) : 0.0); }
    }
    catch (const std::exception& e) { std::cerr << "gpwprobe " << cmd << ": " << e.what() << std::endl; return 1; }
    Usage(); return 2;
}

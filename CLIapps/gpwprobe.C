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
//   gpwprobe mno                   MNO_* / GPW_MNO_*                    the MnO AFM-II campaign run (+ FM arm)
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
//   1x1x1 -> -7.11506  (SiliconGammaConverges)   2x1x1 -> -7.45294  (SiliconMultiKPlumbing, KP-0 re-pin)
//   2x2x2 -> -7.77846  (SR_2x2x2GammaCentred_vs_CP2K)     (2x2x1 has no banked counterpart: timing only)
// SI_XC=becke forces the atom-centred mesh -- the path that exercises SiteStabilizer in the supercell
// setting (the §6a W2b site-adapted angular sets).  Setting cellKind ALONE is a trap: ask for the RECIPE
// (BeckeXCParams), which also makes GPW_BECKE_L / GPW_BECKE_NR live.  The banked anchors are UNIFORM-mesh
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
// THE MnO AFM-II CAMPAIGN RUN (doc/SymmetryUpgradePlan.md §7; doc/Benchmark.md's MnO rows).  Rocksalt MnO
// with the type-II AFM order in the rhombohedral 2-f.u. cell (the FCC cell doubled along [111]).  The CELL
// is the materials entry's named construction; the DISCRIMINATORS are geometry knobs on top of it:
//   MNO_SWAP_SUBLATTICE  put the -m flip on the FIRST Mn (site-exchange equivariance)
//   MNO_SWAP_ORDER       add the (1/2,1/2,1/2) Mn FIRST (position vs atom-index)
//   MNO_SHIFT=f          rigid translation by (f,f,f) fractional (an exact symmetry: everything invariant)
//   MNO_KMESH=n          an n^3 Γ-centred mesh on the magnetic cell (the ordering question needs k)
// The RECIPE knobs: MNO_ORTHO_TOL, MNO_CUTOFF_FACTOR, MNO_ECUT, MNO_SHARED_MU, MNO_MOM_SEED, MNO_REAL,
// MNO_IMPOSE=0/1/2 (free / Shubnikov / grey control), MNO_XC_UNIFORM, MNO_NR, MNO_L, MNO_ALPHA, MNO_KERKER_G0,
// MNO_XC_CUSP, MNO_PULAY, MNO_PULAY_START, MNO_MOM, MNO_MOM_START, MNO_MOM_PENALTY, MNO_MOM_HOLD, MNO_KT,
// GPW_MNO_NMAX, GPW_MNO_VERBOSE; the SCHEDULE: MNO_ANNEAL=kT,kT,... MNO_ACC=... MNO_ANNEAL_PENALTY=...;
// the ARMS: MNO_SKIP_AFM (FM only), MNO_SKIP_FM (AFM only).  Oracle: CP2K MnO AFM-II E=-61.470570 Ha
// (deck IntegrationTests/CP2K/mno_afm2_gpw_sr.inp), Mulliken site moments Mn +/-4.654.
//========================================================================================================
MnOArm RunMnO(int multiplicity, bool afm, const std::string& label)
{
    const double a=8.40;                                       // rocksalt a ~ 4.445 A (a.u.)
    const bool swapSub = std::getenv("MNO_SWAP_SUBLATTICE")!=nullptr;
    const bool swapOrder = std::getenv("MNO_SWAP_ORDER")!=nullptr;
    const double sh = Envd("MNO_SHIFT", 0.0);
    // The named construction (materials.json MnO_AFM2): CubicF re-based by T = FCC doubled along [111].
    auto cellp=std::make_shared<UnitCell>(BravaisCell(Bravais::CubicF, {.a=a}, Matrix3D<int>(0,1,1, 1,0,1, 1,1,0)));
    UnitCell& cell=*cellp;
    if (swapOrder)
    {
        cell.AddAtom(25, {0.5+sh,0.5+sh,0.5+sh}, afm && !swapSub); // -m sublattice, added FIRST
        cell.AddAtom(25, {sh,    sh,    sh    }, afm && swapSub);
    }
    else
    {
        cell.AddAtom(25, {sh,    sh,    sh    }, afm && swapSub);  // Mn sublattice +m (-m when swapped)
        cell.AddAtom(25, {0.5+sh,0.5+sh,0.5+sh}, afm && !swapSub); // Mn sublattice -m (flipped for the AFM arm)
    }
    cell.AddAtom(8,  {0.25+sh,0.25+sh,0.25+sh});
    cell.AddAtom(8,  {0.75+sh,0.75+sh,0.75+sh});
    const int nk=Envi("MNO_KMESH",1);
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));

    MnOArm arm;
    arm.cell=cellp;
    GpwReport report(label, false);

    SolidCalcOptions o;
    o.label=label;
    o.Nelec=26; o.multiplicity=multiplicity;      // 2 x (Mn q7 + O q6); AFM = the two-channel singlet
    o.species={{"Mn",7},{"O",6}};
    o.seed=ChargeDensity::SeedStrategy::IonicSAD;   // Mn2+ d^5 + diffuse O2- -- the basin chooser
    o.ortho=qchem::CholeskyPivoted;                 // cond(S)~7e8: plain Cholesky explodes
    o.orthoTol=Envd("MNO_ORTHO_TOL",1e-4);
    o.cutoffFactor=Envd("MNO_CUTOFF_FACTOR",2.0);
    o.densityEcut =Envd("MNO_ECUT",-1.0);           // <0 = AUTO (cutoffFactor*alpha_max)
    o.spinsShareFermi=Envi("MNO_SHARED_MU",0)!=0;
    o.momFromSeed    =Envi("MNO_MOM_SEED",0)!=0;
    o.forceComplex   =Envi("MNO_REAL",1)==0;        // MNO_REAL=0 -> the all-complex twin
    {   // MNO_IMPOSE: 0 FREE, 1 the SHUBNIKOV group of the declared ordering, 2 the grey erasure control
        const int iv=Envi("MNO_IMPOSE",0);
        o.imposeSymmetry = iv!=0;
        o.greyImposition = iv==2;
    }
    if (std::getenv("MNO_XC_UNIFORM")) o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    if (const char* nr=std::getenv("MNO_NR")) o.xcMesh.nRadial=std::atoi(nr);
    if (const char* ll=std::getenv("MNO_L"))  o.xcMesh.angularDegree=std::atoi(ll);

    SCFParams base;
    base.Verbose=(bool)std::getenv("GPW_MNO_VERBOSE");
    base.StartingRelaxRo=Envd("MNO_ALPHA",0.45);  base.KerkerG0=Envd("MNO_KERKER_G0",1.0);
    base.XCCuspDeficit  =Envd("MNO_XC_CUSP",0.0)!=0.0;
    base.PulayDepth=(int)Envd("MNO_PULAY",0.0);  base.PulayStart=(int)Envd("MNO_PULAY_START",5.0);
    base.UseMOM=Envi("MNO_MOM",1)!=0;  base.MOMStartIter=Envi("MNO_MOM_START",10);
    base.MOMSmearPenalty=Envd("MNO_MOM_PENALTY",0.0);
    base.Guard.HolePersistence=Envi("MNO_MOM_HOLD",3);
    base.SmearingkT=Envd("MNO_KT",5e-3);
    base.NMaxIter=Envi("GPW_MNO_NMAX",80);
    base.MinΔρ=1e-5; base.MinΔE=1e30; base.MinΔFD=1e30; base.MinVirial=1e30; base.MinFD=1e30;
    base.MergeTol=1e-4;

    o.onIteration=[&arm](const SCFIterator::SCFProgress& p)
                  { arm.series.push_back({p.iteration,p.energy,p.dE,p.commutator,p.drho,p.order,p.eb["Eee"]}); };

    auto split=[](const char* env)
    {
        std::vector<std::string> out;
        if (const char* v=std::getenv(env))
            for (std::string t(v), tok; !t.empty(); )
            {
                size_t c=t.find(','); tok=t.substr(0,c);
                if (!tok.empty()) out.push_back(tok);
                if (c==std::string::npos) break;
                t=t.substr(c+1);
            }
        return out;
    };
    auto accType=[](const std::string& n)
    {
        using T=SCFAccelerators::Type;
        if (n=="DIIS")   return T::DIIS;
        if (n=="GDM")    return T::GDM;
        if (n=="Ladder") return T::Ladder;
        if (n=="Null")   return T::Null;
        throw std::runtime_error("MNO_ACC: unknown accelerator \"" + n + "\" (DIIS|GDM|Ladder|Null)");
    };
    const std::vector<std::string> kTs=split("MNO_ANNEAL"), accs=split("MNO_ACC"), lams=split("MNO_ANNEAL_PENALTY");
    if (!(lams.empty() || lams.size()==kTs.size())) throw std::runtime_error("MNO_ANNEAL_PENALTY must parallel MNO_ANNEAL");
    if (!(accs.size()<=1 || kTs.empty() || accs.size()==kTs.size())) throw std::runtime_error("MNO_ACC must parallel MNO_ANNEAL");

    std::vector<SCFStage> schedule;
    if (kTs.empty())
        schedule.push_back({base, accs.empty() ? SCFAccelerators::Type::Ladder : accType(accs[0])});
    else
        for (size_t i=0;i<kTs.size();++i)
        {
            SCFParams p=base;
            p.SmearingkT=std::atof(kTs[i].c_str());
            if (!lams.empty()) p.MOMSmearPenalty=std::atof(lams[i].c_str());
            p.StopOnAccelExhausted = (i+1<kTs.size());
            schedule.push_back({p, accs.empty() ? SCFAccelerators::Type::Ladder
                                  : accType(accs.size()==1 ? accs[0] : accs[i])});
        }

    for (const auto& st : schedule) arm.stageKT.push_back(st.params.SmearingkT);
    arm.calc=std::make_unique<SolidCalculation>(lat, MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR), o, schedule);
    arm.result=arm.calc->Result();
    report::EmitTimings();   // sorted by cost + PEAK RSS, inside the bracket
    return arm;
}

int MnO()
{
    if (std::getenv("MNO_SKIP_AFM"))
    {
        MnOArm F=RunMnO(/*multiplicity*/11, /*afm*/false, "MnO FM Gamma");
        Instrumentation(F, "MnO FM Gamma");
        Check(bool(F.result), "FM converged" + (F.result ? std::string() : ": "+F.result.Error().details));
        if (F.result) CheckNear(F.result->TotalCharge(), 26.0, 1e-6, "FM charge");
        return nFailed;
    }
    MnOArm A=RunMnO(/*multiplicity*/1, /*afm*/true, "MnO AFM-II Gamma");
    Instrumentation(A, "MnO AFM-II Gamma");
    Check(bool(A.result), "AFM-II converged" + (A.result ? std::string() : ": "+A.result.Error().details));
    if (!A.result) return nFailed;
    const ScalarFunction<double>* m=A.result->SpinDensity();
    Check(m!=nullptr, "a multiplicity>=1 run must produce a spin density");
    if (!m) return nFailed;

    // The point probe m(r) at 0.7 bohr off each Mn is a spin DENSITY, not a moment (feedback: integrated
    // observables) -- valid as a COLLAPSE DETECTOR only; the integrated site moment is the m_site column.
    const double a=8.40;
    const rvec3_t off(0.7,0,0), rMn1(0,0,0), rMn2(a,a,a);   // cartesian: A*(1/2,1/2,1/2) = a(1,1,1)
    const double m1=(*m)(rMn1+off), m2=(*m)(rMn2+off);
    std::cout << "[MnO AFM-II] site spin density m(r): Mn1(+seed)="<<m1<<"  Mn2(-seed)="<<m2
              << "  m_stag=½(m1−m2)="<<0.5*(m1-m2)<<"  m_net=m1+m2="<<(m1+m2)
              << (std::abs(m1+m2) > 0.2*std::abs(m1-m2) ? "  ** NOT STAGGERED: the moment sits on ONE sublattice" : "")
              << std::endl;
    CheckNear(A.result->TotalCharge(), 26.0, 1e-6, "AFM-II charge");
    Check(m1 >  0.01, "Mn1 stays in the +m basin the seed chose");
    Check(m2 < -0.01, "Mn2 stays in the -m basin the seed chose");
    Check(std::fabs(m1+m2) <= 0.2*std::fabs(m1), "the two sublattices stagger symmetrically");
    {
        ChargeDensity::fitbasis_t fit(A.calc->Basis().CreateVxcFitBasisSet(A.cell.get(), qcMesh::MeshParams{}));
        auto* pol=dynamic_cast<const ChargeDensity::cSpinResolved_CD*>(&A.result->DensityMatrix());
        Check(pol!=nullptr, "a multiplicity>=1 run must produce a spin-resolved density");
        if (pol)
        {
            auto* fu =dynamic_cast<const ChargeDensity::FourierDensity*>(pol->GetChannel(Spin::Up));
            auto* fdn=dynamic_cast<const ChargeDensity::FourierDensity*>(pol->GetChannel(Spin::Down));
            if (fu && fdn)
            {
                ΔG_Map mu=fu->GetFourierDensity(*fit), md=fdn->GetFourierDensity(*fit);
                const ivec3_t q(1,0,0);
                if (mu.count(q) && md.count(q))
                {
                    const double Mstag=std::abs(dcmplx(mu.at(q))-dcmplx(md.at(q)))*A.cell->GetCellVolume()/2;
                    std::cout << "[MnO AFM-II] |m-tilde(q_AFM)|*Omega/2 = "<<Mstag<<" e- (staggered moment scale)"<<std::endl;
                    Check(Mstag > 1.0, "the converged staggered moment is electrons-scale (d^5 ~ 4-5)");
                }
            }
        }
    }
    if (std::getenv("MNO_SKIP_FM")) return nFailed;
    MnOArm F=RunMnO(/*multiplicity*/11, /*afm*/false, "MnO FM Gamma");
    Instrumentation(F, "MnO FM Gamma");
    Check(bool(F.result), "FM converged" + (F.result ? std::string() : ": "+F.result.Error().details));
    if (!F.result) return nFailed;
    CheckNear(F.result->TotalCharge(), 26.0, 1e-6, "FM charge");
    const double Eafm=A.result->Energy(), Efm=F.result->Energy();
    std::cout << "[MnO ordering] E_AFM="<<std::setprecision(10)<<Eafm<<"  E_FM="<<Efm<<"  dE="<<(Efm-Eafm)*1000<<" mHa"<<std::endl;
    Check(Eafm < Efm, "AFM-II is the LSDA ground-state ordering");
    return nFailed;
}

void Usage()
{
    std::cout << "gpwprobe ladder | ksweep | naf-smear | becke-ladder {si|naf|mn|al} | mno\n"
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
        if (cmd=="mno")          return MnO();
    }
    catch (const std::exception& e) { std::cerr << "gpwprobe " << cmd << ": " << e.what() << std::endl; return 1; }
    Usage(); return 2;
}

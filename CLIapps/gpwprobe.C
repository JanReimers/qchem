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
//   gpwprobe nio                   NIO_* / GPW_NIO_*                    the SAME arm on NiO (the ACBN0-vs-hp.x cell)
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
// THE ROCKSALT AFM-II ARM'S MATERIAL SPEC.  MnO and NiO are the SAME run -- the same cell construction, the
// same ionic-SAD basin choice, the same knobs -- one lattice constant and one transition metal apart, so
// they are one arm with two specs rather than two copies.  Each spec owns its ENV PREFIX: MnO's knobs stay
// MNO_* (every banked recipe in doc/Benchmark.md is written that way and must keep working verbatim) and
// NiO's are NIO_*, with no fallback between them -- a run says in its own knob names which material it drove.
//========================================================================================================
struct TmoSpec
{
    std::string prefix;      //!< env-var prefix: "MNO" | "NIO"
    std::string name;        //!< "MnO" | "NiO"; the cell mirrors materials.json's <name>_AFM2 entry
    std::string tm;          //!< the transition metal's symbol, "Mn" | "Ni"
    int    Z=0;              //!< the transition metal's atomic number
    int    qTM=0, qO=6;      //!< GTH valence electrons (Mn q7 / Ni q10; O q6)
    double a=0.0;            //!< the re-based cell's cubic lattice constant (a.u.)
    int    Nelec=0;          //!< 2*(qTM+qO)
    int    multFM=0;         //!< the FM arm's multiplicity (d^5 -> 11, d^8 -> 5)
    ChargeDensity::SeedStrategy seed=ChargeDensity::SeedStrategy::IonicSAD;  //!< see the two specs below
    //! Env lookup under this spec's prefix: Env("KMESH") reads MNO_KMESH or NIO_KMESH.
    const char* Env(const char* knob) const { return std::getenv((prefix+"_"+knob).c_str()); }
    double Envd(const char* knob, double d) const { const char* s=Env(knob); return s ? std::atof(s) : d; }
    int    Envi(const char* knob, int    d) const { const char* s=Env(knob); return s ? std::atoi(s) : d; }
};
const TmoSpec MnOSpec{"MNO","MnO","Mn",25, 7,6, 8.40, 26, 11, ChargeDensity::SeedStrategy::IonicSAD};
// NiO's a = 7.88 is QE's hp.x benchmark cell, which is where U(Ni 3d) = 5.27 eV was banked: the estimator
// must be compared with the oracle on the oracle's geometry (doc/OpenWork.md step 5, increment 3).
// ⚠ **ITS SEED WAS SAD FOR A DAY, ON A CONCLUSION THAT WAS WRONG** (retracted 2026-09-23).  The claim was
// that Ni2+ cannot be generated by our pseudo-atom -- high-spin d^8 being minority-d^3 in a five-fold shell,
// with the atom occupying whole angular irreps -- on evidence of a NON-AUFBAU run whose density ran away to
// <r> = 2.91 bohr for a CATION.  That run predated the iteration-cap fix in the SAME session: the pseudo-atom
// was simply still descending at the default 20 iterations.  At --nmax 120 Ni2+ converges cleanly, charge
// 8.001, moment 2.000, <r> = 1.008 bohr -- between Mn3+ (1.078) and Ni3+ (0.938), exactly where a d^8 cation
// belongs.  So IonicSAD is available and NiO uses it, like MnO.
// ★ What IS true, and much milder: a PARTIALLY FILLED MINORITY shell under a filled majority is a LONG
// DESCENT, because two nearly-degenerate orbitals straddle the occupied/empty boundary (eps gap 2e-6 Ha at
// iteration 60).  It converges; it just needs the iterations.  A partially filled MAJORITY over an empty
// minority does not have the problem at all -- Mn3+ (d^4), Mn4+ (d^3), Ni3+ (d^7) and Co3+ (d^6) all converge
// by 60.  ⇒ the lesson is the cap one (a converged=false past the 3d row is the CAP until proven otherwise),
// not a structural limit of the occupation recipe.
const TmoSpec NiOSpec{"NIO","NiO","Ni",28,10,6, 7.88, 32,  5, ChargeDensity::SeedStrategy::IonicSAD};

//========================================================================================================
// THE ROCKSALT AFM-II CAMPAIGN RUN (doc/SymmetryUpgradePlan.md §7; doc/Benchmark.md's MnO rows).  A rocksalt
// transition-metal monoxide with the type-II AFM order in the rhombohedral 2-f.u. cell (the FCC cell doubled
// along [111]).  TWO MATERIALS, ONE BODY: `gpwprobe mno` is the MnO campaign, `gpwprobe nio` the same run on
// NiO -- the cell the hp.x oracle used for U(Ni 3d) = 5.27 eV (doc/OpenWork.md step 5, increment 3).  Knobs
// below are written <P>_ for the spec's prefix: MNO_ for mno, NIO_ for nio, no fallback between them.  The
// CELL is the materials entry's named construction; the DISCRIMINATORS are geometry knobs on top of it:
//   <P>_SWAP_SUBLATTICE  put the -m flip on the FIRST TM site (site-exchange equivariance)
//   <P>_SWAP_ORDER       add the (1/2,1/2,1/2) TM site FIRST (position vs atom-index)
//   <P>_SHIFT=f          rigid translation by (f,f,f) fractional (an exact symmetry: everything invariant)
//   <P>_KMESH=n          an n^3 Γ-centred mesh on the magnetic cell (the ordering question needs k)
// The RECIPE knobs: <P>_ORTHO_TOL, <P>_CUTOFF_FACTOR, <P>_ECUT, <P>_SHARED_MU, <P>_MOM_SEED, <P>_REAL,
// <P>_IMPOSE=0/1/2 (free / Shubnikov / grey control), <P>_XC_UNIFORM, <P>_XC_ECUT=Ha, <P>_VET=1 (the pin-22 vet-stage basis trim), <P>_NR, <P>_L, <P>_ALPHA, <P>_KERKER_G0,
// <P>_XC_CUSP, <P>_PULAY, <P>_PULAY_START, <P>_MOM, <P>_MOM_START, <P>_MOM_PENALTY, <P>_MOM_HOLD, <P>_KT,
// GPW_<P>_NMAX, GPW_<P>_VERBOSE, <P>_CHI0=nq (chi0 on an nq^3 q-mesh, LinearResponsePlan R0), <P>_CHI=1 (the
// self-consistent chi, chi0 and U at q=0, R2 -- needs <P>_REAL=0), <P>_U=eV (DFT+U on both TM d, programme step 5) + <P>_U_IRREP=a,b,c (eV per
// site-irrep slot, increment 2: a1g<t2g, e_g<e_g, e_g<t2g under D_3d) + <P>_ACBN0=1 (print the ACBN0 (U,J)
// estimate from the converged orbitals; <P>_ACBN0=n>1 runs the paper's OUTER LOOP for up to n steps, re-converging
// on the same Hamiltonian, <P>_ACBN0_TOL=eV; increment 3) + <P>_U_RADIAL=every|atomic|ortho|orthofull (the
// manifold's radial and projector: CP2K's every shell, the pseudo-atom 3d, ortho-atomic among the TM, or QE's
// full ortho-atomic set with O 2s/2p + TM 4s spectators at U=0), <P>_EPS=tol +
// <P>_MEASURE=maxdd|mixer (CP2K's EPS_SCF measure max|dD_ij| between successive D_out, or the mixer's own
// residual -- doc/Benchmark.md rule 3f: an iteration count is comparable only on the same measure); the SCHEDULE: <P>_ANNEAL=kT,kT,... <P>_ACC=... <P>_ANNEAL_PENALTY=...;
// the ARMS: <P>_SKIP_AFM (FM only), <P>_SKIP_FM (AFM only); the STATE (CK-1): <P>_SAVE=path (write the state
// after every stage), <P>_RESTART=path (start from one, converging with the schedule's final stage); the FM arm
// reads/writes path.fm.  Oracle: CP2K MnO AFM-II E=-61.470570 Ha
// (deck IntegrationTests/CP2K/mno_afm2_gpw_sr.inp), Mulliken site moments Mn +/-4.654.  NiO's oracle is hp.x,
// not CP2K: there is no CP2K deck for it and no banked total.
//========================================================================================================
MnOArm RunTMO(const TmoSpec& S, int multiplicity, bool afm, const std::string& label)
{
    const double a=S.a;
    const bool swapSub = S.Env("SWAP_SUBLATTICE")!=nullptr;
    const bool swapOrder = S.Env("SWAP_ORDER")!=nullptr;
    const double sh = S.Envd("SHIFT", 0.0);
    // The named construction (materials.json <name>_AFM2): CubicF re-based by T = FCC doubled along [111].
    auto cellp=std::make_shared<UnitCell>(BravaisCell(Bravais::CubicF, {.a=a}, Matrix3D<int>(0,1,1, 1,0,1, 1,1,0)));
    UnitCell& cell=*cellp;
    if (swapOrder)
    {
        cell.AddAtom(S.Z, {0.5+sh,0.5+sh,0.5+sh}, afm && !swapSub); // -m sublattice, added FIRST
        cell.AddAtom(S.Z, {sh,    sh,    sh    }, afm && swapSub);
    }
    else
    {
        cell.AddAtom(S.Z, {sh,    sh,    sh    }, afm && swapSub);  // TM sublattice +m (-m when swapped)
        cell.AddAtom(S.Z, {0.5+sh,0.5+sh,0.5+sh}, afm && !swapSub); // TM sublattice -m (flipped for the AFM arm)
    }
    cell.AddAtom(8,  {0.25+sh,0.25+sh,0.25+sh});
    cell.AddAtom(8,  {0.75+sh,0.75+sh,0.75+sh});
    const int nk=S.Envi("KMESH",1);
    Lattice_3D lat(cell, ivec3_t(nk,nk,nk));

    MnOArm arm;
    arm.cell=cellp;
    GpwReport report(label, false);

    SolidCalcOptions o;
    o.label=label;
    o.Nelec=S.Nelec; o.multiplicity=multiplicity;   // 2 x (TM + O q6); AFM = the two-channel singlet
    o.species={{S.tm,S.qTM},{"O",S.qO}};
    o.seed=S.seed;                                  // IonicSAD (TM2+ d^n + diffuse O2-) for MnO; see NiOSpec
    o.ortho=qchem::CholeskyPivoted;                 // cond(S)~7e8: plain Cholesky explodes
    o.orthoTol=S.Envd("ORTHO_TOL",1e-4);
    o.cutoffFactor=S.Envd("CUTOFF_FACTOR",2.0);
    o.densityEcut =S.Envd("ECUT",-1.0);           // <0 = AUTO (cutoffFactor*alpha_max)
    o.spinsShareFermi=S.Envi("SHARED_MU",0)!=0;
    o.momFromSeed    =S.Envi("MOM_SEED",0)!=0;
    o.forceComplex   =S.Envi("REAL",1)==0;        // <P>_REAL=0 -> the all-complex twin
    {   // <P>_IMPOSE: 0 FREE, 1 the SHUBNIKOV group of the declared ordering, 2 the grey erasure control
        const int iv=S.Envi("IMPOSE",0);
        o.imposeSymmetry = iv!=0;
        o.greyImposition = iv==2;
    }
    if (S.Env("XC_UNIFORM")) o.xcMesh.cellKind=qcMesh::UnitCellKind::Uniform;
    if (const char* ec=S.Env("XC_ECUT")) o.xcMesh.eCut=std::atof(ec);   // the uniform XC mesh's cutoff (Ha); 0 = the manual nUniform
    // <P>_NR / <P>_L are BECKE-mesh knobs.  Under the default Auto choice the facade rebuilds the Becke
    // parameters from scratch, so setting nRadial/angularDegree on an Auto mesh was silently ignored
    // (found 2026-09-27, the NiO Becke-convergence arm): an explicit NR or L now PINS the Becke mesh with them.
    if (!S.Env("XC_UNIFORM") && (S.Env("NR") || S.Env("L")))
    {
        o.xcMesh=qcMesh::BeckeXCParams(S.Envi("NR",-1), -1.0, S.Envi("L",-1));
        o.xcMesh.cellKind=qcMesh::UnitCellKind::Becke;
    }
    else
    {
        if (const char* nr=S.Env("NR")) o.xcMesh.nRadial=std::atoi(nr);
        if (const char* ll=S.Env("L"))  o.xcMesh.angularDegree=std::atoi(ll);
    }
    // <P>_ACBN0=1 carries the TM d manifolds even at U=0 and, after the arm converges, prints the ACBN0
    // estimate (U-bar, J-bar, U_eff) from its orbitals -- one step of the paper's outer loop; iterate by
    // hand with <P>_U=<U_eff of the previous run> (increment 3).
    const bool acbn0 = S.Envi("ACBN0",0)!=0;
    if (const double U=S.Envd("U",0.0); U>0.0 || acbn0 || S.Envi("CHI0",0)>0 || S.Envi("CHI",0)>0)   // CHI0/CHI: the channels ARE the manifolds
    {
        // <P>_U_RADIAL: every (CP2K's every-shell manifold, default) | atomic (ONE contracted pseudo-atom 3d, QE's
        // `atomic`) | ortho (the two TM 3d sets Löwdin-orthogonalised against each other) | orthofull (QE's
        // ortho-atomic SET: + O 2s, O 2p and TM 4s as U=0 spectators, sites 2,3 = O) -- increment 3 slices C/D.
        const std::string rad = S.Env("U_RADIAL") ? S.Env("U_RADIAL") : "every";
        if      (rad=="every")  o.hubbard={HubbardU(0,2,U), HubbardU(1,2,U)};                      // sites 0,1 = the TM
        else if (rad=="atomic") o.hubbard={HubbardU_Atomic(0,2,U), HubbardU_Atomic(1,2,U)};
        else if (rad=="ortho")  o.hubbard={HubbardU_OrthoAtomic(0,2,U), HubbardU_OrthoAtomic(1,2,U)};
        else if (rad=="orthofull") o.hubbard={HubbardU_OrthoAtomic(0,2,U), HubbardU_OrthoAtomic(1,2,U),
                                              HubbardU_OrthoAtomic(0,0,0.0), HubbardU_OrthoAtomic(1,0,0.0),      // TM 4s
                                              HubbardU_OrthoAtomic(2,0,0.0), HubbardU_OrthoAtomic(3,0,0.0),      // O 2s
                                              HubbardU_OrthoAtomic(2,1,0.0), HubbardU_OrthoAtomic(3,1,0.0)};     // O 2p
        else throw std::runtime_error(S.prefix+"_U_RADIAL: expected every|atomic|ortho|orthofull, got '"+rad+"'");
        // <P>_U_IRREP=a,b,c (eV): one U per site-irrep slot, in the order of the term's "[+U] site .. U slots"
        // table -- for the AFM-II Mn under D_3d < O_h that is [0] a1g<t2g, [1] e_g<e_g, [2] e_g<t2g (increment 2).
        // The count must match the slot count or the term throws with the table in the message.
        if (const char* v=S.Env("U_IRREP"))
            for (std::string t(v), tok; !t.empty(); )
            {
                size_t c=t.find(','); tok=t.substr(0,c);
                if (!tok.empty()) for (auto& M : o.hubbard) M.Uirrep.push_back(std::stod(tok)/27.211386245988);
                if (c==std::string::npos) break;
                t=t.substr(c+1);
            }
    }

    SCFParams base;
    base.Verbose=std::getenv(("GPW_"+S.prefix+"_VERBOSE").c_str())!=nullptr;
    base.StartingRelaxRo=S.Envd("ALPHA",0.45);  base.KerkerG0=S.Envd("KERKER_G0",1.0);
    base.XCCuspDeficit  =S.Envd("XC_CUSP",0.0)!=0.0;
    base.PulayDepth=(int)S.Envd("PULAY",0.0);  base.PulayStart=(int)S.Envd("PULAY_START",5.0);
    base.UseMOM=S.Envi("MOM",1)!=0;  base.MOMStartIter=S.Envi("MOM_START",10);
    base.MOMSmearPenalty=S.Envd("MOM_PENALTY",0.0);
    base.Guard.HolePersistence=S.Envi("MOM_HOLD",3);
    base.SmearingkT=S.Envd("KT",5e-3);
    // 200, not 80: the old default was tuned on MnO, which converges in ~18 iterations, and NiO at U=3 needs
    // 114 -- so 80 reported "NOT converged" on a run that was descending perfectly well (2026-09-23, the THIRD
    // time an iteration cap produced a wrong conclusion this week).  A cap only costs anything when it is HIT,
    // so buy margin: it changes no run that was already converging inside 80.
    base.NMaxIter=[&]{ const char* v=std::getenv(("GPW_"+S.prefix+"_NMAX").c_str()); return v?std::atoi(v):200; }();
    base.MinΔρ=S.Envd("EPS",1e-5); base.MinΔE=1e30; base.MinΔFD=1e30; base.MinVirial=1e30; base.MinFD=1e30;
    if (const char* m=S.Env("MEASURE"))
    {
        const std::string ms(m);
        if      (ms=="maxdd") base.Δρmeasure=SCFParams::Measure::MaxΔD;
        else if (ms=="mixer") base.Δρmeasure=SCFParams::Measure::MixerResidual;
        else throw std::runtime_error(S.prefix+"_MEASURE: expected maxdd|mixer, got '"+ms+"'");
    }
    base.MergeTol=1e-4;

    o.onIteration=[&arm](const SCFIterator::SCFProgress& p)
                  { arm.series.push_back({p.iteration,p.energy,p.dE,p.commutator,p.drho,p.order,p.eb["Eee"]}); };

    auto split=[&S](const char* knob)
    {
        std::vector<std::string> out;
        if (const char* v=S.Env(knob))
            for (std::string t(v), tok; !t.empty(); )
            {
                size_t c=t.find(','); tok=t.substr(0,c);
                if (!tok.empty()) out.push_back(tok);
                if (c==std::string::npos) break;
                t=t.substr(c+1);
            }
        return out;
    };
    auto accType=[&S](const std::string& n)
    {
        using T=SCFAccelerators::Type;
        if (n=="DIIS")   return T::DIIS;
        if (n=="GDM")    return T::GDM;
        if (n=="Ladder") return T::Ladder;
        if (n=="Null")   return T::Null;
        throw std::runtime_error(S.prefix+"_ACC: unknown accelerator \"" + n + "\" (DIIS|GDM|Ladder|Null)");
    };
    const std::vector<std::string> kTs=split("ANNEAL"), accs=split("ACC"), lams=split("ANNEAL_PENALTY");
    if (!(lams.empty() || lams.size()==kTs.size())) throw std::runtime_error(S.prefix+"_ANNEAL_PENALTY must parallel "+S.prefix+"_ANNEAL");
    if (!(accs.size()<=1 || kTs.empty() || accs.size()==kTs.size())) throw std::runtime_error(S.prefix+"_ACC must parallel "+S.prefix+"_ANNEAL");

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
    // ⚠ A SPHERICAL-d arm needs the VA/SPH span (user, 2026-09-21): the SR Mn block carries only TWO s exponents
    // because its s span comes from the CARTESIAN d shells' x^2+y^2+z^2 contaminants -- under GPW_SPHERICAL=1
    // those are gone.  A +U manifold needs the spherical view, so an unset GPW_BASIS_SPAN defaults to VA here
    // (the oracle gate's span); an explicit GPW_BASIS_SPAN=sr with GPW_SPHERICAL is refused.  Ni is in sr for
    // no arm -- its block exists only in va/sph (BasisSetData/valence_lowq_va.bsd's Ni note) -- so the same
    // default and the same refusal are what NiO needs too.
    if (std::getenv("GPW_SPHERICAL"))
    {
        const char* span=std::getenv("GPW_BASIS_SPAN");
        if (!span) { setenv("GPW_BASIS_SPAN", "va", 1); std::cout << "[" << S.name << "] spherical d => GPW_BASIS_SPAN defaulted to va (the SR " << S.tm << " block has no s span without the Cartesian d contaminants)" << std::endl; }
        else if (std::string(span)=="sr") throw std::runtime_error("gpwprobe: GPW_BASIS_SPAN=sr with GPW_SPHERICAL -- the SR "+S.tm+" block's s span lives in the Cartesian d contaminants; use va or sph");
    }
    // <P>_VET=1: the VET-STAGE trim (doc/Pins.md pin 22) -- near-dependent diffuse shells removed ONCE, per
    // element, on the full k-mesh, BEFORE the run is built, at the run's own orthoTol -- so the ortho step has
    // nothing left to drop in any k-block.  The ortho-time per-k drop it replaces made the basis differ from k
    // to k and between equivalent sites (NiO AFM-II, 2026-09-27).
    std::shared_ptr<const BasisSet::Real_BS> mol;
    if (S.Envi("VET",0)!=0)
    {
        auto make=[&cell](const BasisSet::Gaussian::ShellTrim& t){ return MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR, t); };
        mol=BasisSet::Lattice::VetStageTrim(lat, make, {.images=o.images, .kShift=o.kShift}, o.orthoTol).mol;
    }
    else if (const char* tr=S.Env("TRIM"))
    {   // <P>_TRIM=Z:l:alpha[,Z:l:alpha...] -- a STATED trim, built ONCE with no vet loop in the process (the A/B
        // for anything the vet's trial builds could leave behind)
        BasisSet::Gaussian::ShellTrim t;
        for (std::string rest(tr), tok; !rest.empty(); )
        {
            size_t c=rest.find(','); tok=rest.substr(0,c);
            int Z=0, l=0; double a=0;
            if (std::sscanf(tok.c_str(), "%d:%d:%lf", &Z, &l, &a)!=3) throw std::runtime_error(S.prefix+"_TRIM: expected Z:l:alpha, got '"+tok+"'");
            t.shells.push_back({Z,l,rvec_t(1,a)});
            if (c==std::string::npos) break;
            rest=rest.substr(c+1);
        }
        std::cout << "[basis trim] STATED (no vet loop): "; t.Write(std::cout); std::cout << std::endl;
        mol=MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR, t);
    }
    else mol=MakeBasisLowQ(cell, BasisSetData::VALENCE_LOWQ_SR);
    // <P>_SAVE=path: write the state after every stage (CK-1 saveStateTo -- the last write is the final stage).
    // <P>_RESTART=path: start from a saved state instead of the seed.  Restart runs ONE stage, so it takes the
    // schedule's FINAL stage (params + accelerator): the earlier stages of an annealed recipe exist only to deliver
    // a good starting density, which is what the file already holds (CK-1 residual c: no schedule overload yet).
    // A refusal (foreign file, a REFUSE-class fingerprint difference) THROWS: a probe never silently re-seeds.
    // The FM arm's file is the path + ".fm", so a run of both arms cannot overwrite the AFM state with the FM one.
    const std::string armTag = afm ? "" : ".fm";
    if (const char* sv=S.Env("SAVE")) o.saveStateTo=std::string(sv)+armTag;
    if (const char* rs0=S.Env("RESTART"))
    {
        const std::string rs=std::string(rs0)+armTag;
        SolidCalcOptions ro=o;
        ro.accelerator=schedule.back().accelerator;
        auto R=SolidCalculation::Restart(rs, lat, mol, ro, schedule.back().params);
        if (!R) throw std::runtime_error(S.prefix+"_RESTART="+rs+" refused: "+R.Error().details);
        arm.calc=R.TakeValue();
    }
    else arm.calc=std::make_unique<SolidCalculation>(lat, mol, o, schedule);
    arm.result=arm.calc->Result();
    if (acbn0)
    {
        std::cout << "[" << S.name << " " << o.label << "] ACBN0 from the last iterate"
                  << (arm.result ? " (CONVERGED):" : " (NOT converged -- a diagnostic, not a U):") << std::endl;
        const int nOuter=S.Envi("ACBN0",0);
        if (nOuter<=1) arm.calc->EstimateHubbardU();                  // one-shot
        else
        {   // the paper's outer loop, on the same Hamiltonian, re-converging with the schedule's FINAL stage
            SolidCalculation::HubbardLoop lp; lp.maxOuter=size_t(nOuter); lp.tolU_eV=S.Envd("ACBN0_TOL",1e-3);
            auto R=arm.calc->ConvergeHubbardU(schedule.back().params, lp);
            std::cout << "[" << S.name << " " << o.label << "] ACBN0 loop: " << R.outer << " outer steps, "
                      << (R.converged ? "U CONVERGED" : "U NOT converged") << " (tol " << lp.tolU_eV << " eV), last SCF "
                      << (R.scfConverged ? "converged" : "NOT converged") << ";  U_eff trajectory (site 0, eV):";
            for (const auto& u : R.U_eV) std::cout << " " << u[0];
            std::cout << std::endl;
            arm.result=arm.calc->Result();                             // the final-U SCF is now the arm's answer
        }
    }
    // <P>_CHI0=n: chi0 over the Hubbard manifolds on an n^3 q-mesh (doc/LinearResponsePlan.md stage R0) -- the
    // independent-particle response hp.x prints at its first iteration.  Needs a FULL k-mesh (no <P>_IMPOSE)
    // commensurate with the q-mesh, and the channels listed as +U manifolds (at U=0 to probe without +U).
    if (const int nq=S.Envi("CHI0",0); nq>0)
    {
        std::cout << "[" << S.name << " " << o.label << "] chi0 from the last iterate"
                  << (arm.result ? " (CONVERGED):" : " (NOT converged -- a diagnostic only):") << std::endl;
        auto chi=arm.calc->IndependentResponse(ivec3_t(nq,nq,nq));   // reports itself
        (void)chi;
    }
    // <P>_CHI=1: the SELF-CONSISTENT response at q = 0 over the Hubbard manifolds (doc/LinearResponsePlan.md
    // stage R2): chi0, chi and U = (chi0^-1 - chi^-1)_II of THIS cell (the LR-cDFT U of a supercell of this
    // size).  Needs the complex ansatz (<P>_REAL=0) and a full k-mesh; at <P>_U=0 this is the U_0 of the U=0
    // ground state (+U answers zero); an unfrozen U != 0 is refused.
    if (S.Envi("CHI",0)>0)
    {
        std::cout << "[" << S.name << " " << o.label << "] self-consistent chi (q=0) from the last iterate"
                  << (arm.result ? " (CONVERGED):" : " (NOT converged -- a diagnostic only):") << std::endl;
        auto chi=arm.calc->HubbardLinearResponse();   // reports itself
        (void)chi;
    }
    report::EmitTimings();   // sorted by cost + PEAK RSS, inside the bracket
    return arm;
}

//! The arm's DRIVER: the AFM-II run with its magnetic diagnostics, then (unless skipped) the FM run and the
//! ordering comparison.  One body for both materials -- the checks are statements about a rocksalt AFM-II
//! antiferromagnet, not about Mn.
int TMO(const TmoSpec& S)
{
    const std::string afmLabel=S.name+" AFM-II Gamma", fmLabel=S.name+" FM Gamma";
    if (S.Env("SKIP_AFM"))
    {
        MnOArm F=RunTMO(S, S.multFM, /*afm*/false, fmLabel);
        Instrumentation(F, fmLabel);
        Check(bool(F.result), "FM converged" + (F.result ? std::string() : ": "+F.result.Error().details));
        if (F.result) CheckNear(F.result->TotalCharge(), double(S.Nelec), 1e-6, "FM charge");
        return nFailed;
    }
    MnOArm A=RunTMO(S, /*multiplicity*/1, /*afm*/true, afmLabel);
    Instrumentation(A, afmLabel);
    Check(bool(A.result), "AFM-II converged" + (A.result ? std::string() : ": "+A.result.Error().details));
    if (!A.result) return nFailed;
    const ScalarFunction<double>* m=A.result->SpinDensity();
    Check(m!=nullptr, "a multiplicity>=1 run must produce a spin density");
    if (!m) return nFailed;

    // The point probe m(r) at 0.7 bohr off each TM site is a spin DENSITY, not a moment (feedback: integrated
    // observables) -- valid as a COLLAPSE DETECTOR only; the integrated site moment is the m_site column.
    const double a=S.a;
    const rvec3_t off(0.7,0,0), r1(0,0,0), r2(a,a,a);       // cartesian: A*(1/2,1/2,1/2) = a(1,1,1)
    const double m1=(*m)(r1+off), m2=(*m)(r2+off);
    std::cout << "["<<S.name<<" AFM-II] site spin density m(r): "<<S.tm<<"1(+seed)="<<m1<<"  "<<S.tm<<"2(-seed)="<<m2
              << "  m_stag=½(m1−m2)="<<0.5*(m1-m2)<<"  m_net=m1+m2="<<(m1+m2)
              << (std::abs(m1+m2) > 0.2*std::abs(m1-m2) ? "  ** NOT STAGGERED: the moment sits on ONE sublattice" : "")
              << std::endl;
    CheckNear(A.result->TotalCharge(), double(S.Nelec), 1e-6, "AFM-II charge");
    Check(m1 >  0.01, S.tm+"1 stays in the +m basin the seed chose");
    Check(m2 < -0.01, S.tm+"2 stays in the -m basin the seed chose");
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
                    std::cout << "["<<S.name<<" AFM-II] |m-tilde(q_AFM)|*Omega/2 = "<<Mstag<<" e- (staggered moment scale)"<<std::endl;
                    // The scale is the TM's unpaired-electron count: MnO d^5 ~ 4-5, NiO d^8 ~ 2.  One electron
                    // is the floor BOTH clear when the order survives, and neither clears when it collapses.
                    Check(Mstag > 1.0, "the converged staggered moment is electrons-scale");
                }
            }
        }
    }
    if (S.Env("SKIP_FM")) return nFailed;
    MnOArm F=RunTMO(S, S.multFM, /*afm*/false, fmLabel);
    Instrumentation(F, fmLabel);
    Check(bool(F.result), "FM converged" + (F.result ? std::string() : ": "+F.result.Error().details));
    if (!F.result) return nFailed;
    CheckNear(F.result->TotalCharge(), double(S.Nelec), 1e-6, "FM charge");
    const double Eafm=A.result->Energy(), Efm=F.result->Energy();
    std::cout << "["<<S.name<<" ordering] E_AFM="<<std::setprecision(10)<<Eafm<<"  E_FM="<<Efm<<"  dE="<<(Efm-Eafm)*1000<<" mHa"<<std::endl;
    Check(Eafm < Efm, "AFM-II is the LSDA ground-state ordering");
    return nFailed;
}

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
    std::cout << "gpwprobe ladder | ksweep | naf-smear | becke-ladder {si|naf|mn|al} | mno | nio | gate1 <material> [U_eV]\n"
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
        if (cmd=="mno")          return TMO(MnOSpec);
        if (cmd=="nio")          return TMO(NiOSpec);
        if (cmd=="gate1")        { if (argc<3) { Usage(); return 2; }
                                    return Gate1(argv[2], argc>3 ? std::atof(argv[3]) : 0.0); }
    }
    catch (const std::exception& e) { std::cerr << "gpwprobe " << cmd << ": " << e.what() << std::endl; return 1; }
    Usage(); return 2;
}

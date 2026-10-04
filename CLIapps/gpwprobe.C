// File: CLIapps/gpwprobe.C  THE GPW CAMPAIGN INSTRUMENTS -- the hand-run ladders, sweeps and probes that used to
// sit in IntegrationTests/GPW_SCF_UT.C as DISABLED_ tests (doc/TestSuitePlan.md §8, verdict P, 2026-09-15).
//
// ONE instrument survives here (D-ENV step 6d.5, 2026-10-04):
//
//   gpwprobe becke-ladder SYSTEM   SYSTEM = si | naf | mn | al          the V2.6 Becke (nR, degree) ladder: converge one system, freeze its density, and
//                                                                       score a ladder of Becke rules -- an instrument that READS, not a run to reproduce
//
// Every other probe that lived here (ladder, ksweep, naf-smear, gate1, mno, nio) is now an input DECK (rundeck <deck> [--set key=value]...) or a hard-coded
// integration test; their environment knobs are gone.  The table below says where each went.  An instrument is something a human runs with knobs and
// READS; a test is something ctest runs and JUDGES -- and the judging half moved to ITMain.
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

//========================================================================================================
// RETIRED 2026-10-04 (D-ENV step 6d.5), each as a DECK (rundeck <deck> [--set key=value]...) or a hard-coded test -- the env knobs went with them:
//   ladder     (SI_LADDER, SI_XC)      -> decks/ladder/Si_{2x1x1,2x2x1,2x2x2}_gamma.json (supercells are materials.json entries); the folding GATE is
//                                         ITMain GPW_Si.Γ_Imp_eqK211.  SI_XC=becke is `--set solid.xcMesh.cellKind=Becke`.
//   ksweep     (GPW_KSHIFT)            -> decks/probes/Si_single_k.json  --set solid.kShift=[s,s,s]
//   naf-smear  (NAF_ECUT, NAFGDM_*)    -> decks/probes/NaF_gdm_smear.json
//   gate1      (GATE1_*)               -> decks/probes/{LiMn2O4_spinel,MnO2_lambda_spinel}_gate1.json  --set solid.hubbard.N.U_eV=...
//   mno / nio  (MNO_*, NIO_*)          -> decks/MnO_AFM2_free.json, decks/NiO_AFM2_acbn0.json   (6d.3)
// What stays is the one INSTRUMENT that is not a run: becke-ladder scores quadrature rules on a frozen density (no knobs, nothing to reproduce from a deck).
//========================================================================================================

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

void Usage()
{
    std::cout << "gpwprobe becke-ladder {si|naf|mn|al}\n  (every other probe is retired: they are decks -- see the table at the top of CLIapps/gpwprobe.C)\n";
}
} // anonymous

int main(int argc, char** argv)
{
    if (argc<2) { Usage(); return 2; }
    const std::string cmd=argv[1];
    try
    {
        if (cmd=="becke-ladder") { if (argc<3) { Usage(); return 2; } return BeckeRecipeLadder(argv[2]); }
        if (cmd=="ladder" || cmd=="ksweep" || cmd=="naf-smear" || cmd=="gate1" || cmd=="mno" || cmd=="nio")
        { std::cerr << "gpwprobe " << cmd << " is RETIRED (D-ENV 6d): it is a DECK now -- see the table at the top of CLIapps/gpwprobe.C and doc/Records/EnvKnobInventory.md §10" << std::endl; return 2; }
    }
    catch (const std::exception& e) { std::cerr << "gpwprobe " << cmd << ": " << e.what() << std::endl; return 1; }
    Usage(); return 2;
}

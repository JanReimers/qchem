// File: Calculation/Imp/SolidSetup.C  The periodic run setup: every decision made ONCE, before a geometry point is built.
module;
#include <algorithm>   // std::max (alpha_pp)
#include <cassert>
#include <cmath>       // std::pow (the strain allowance's cost)
#include <iomanip>
#include <sstream>
#include <iostream>    // the stated-trim line
#include <memory>
#include <optional>    // the what-if's recipe
#include <stdexcept>
#include <string>
#include <tuple>       // std::tie (decision 2)
#include <utility>
#include <vector>
module qchem.SolidSetup;

import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;   // ToleranceBinding (the recipe binds the run's pair-loop tolerances)
import qchem.BasisSet.Gaussian.Lattice.RasterLadder;   // ExponentAtCutoff, kResolutionEps: alpha_max from the basis's Resolution
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH local PP -> alpha_pp
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.Mesh.XCPolicy;                    // XCMeshSharpness / ResolveXCMesh (the grid decision)
import qchem.SeedPolicy;   // ResolveSeed: the ONE Auto seed rule (shared with MoleculeSetup)
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView (GaussianBasisRecipe::angular)
import qchem.Policy;                           // policy::Refuse: each refusal files its row where it fires
import qchem.Reporting;                         // report::Scope: the vet-stage trim files under basisTrim
import qchem.XCRhoPolicy;                      // ReportSolidXCRho: the XC rho-axis decision
import qchem.BasisCoverage;                    // RequireBasisCovers: a basis lacking an element is refused before the build

namespace qchem
{

//! The seed library key for this run: the set's functional and, for each element of the structure, the
//! GTH variant q its Hamiltonian uses -- so the seed density is the SAME pseudo-atom as the terms (C4 step 3).
static ChargeDensity::SeedKey SeedKeyOf(const Pseudopotential::PseudopotentialSet& pp, const Structure& st)
{
    ChargeDensity::SeedKey key(pp.Functional());
    for (const auto& [element, q] : pp.Species(st)) key.q[(int)thePeriodicTable().GetZ(element)] = q;
    return key;
}

//---------------------------------------------------------------------------------------------------
//  The run's SHARPNESS -- an above-SCFIterator decision input.
//
//  Both sources come off ABSTRACT capability faces via the sanctioned abstract->abstract cross-cast, so
//  nothing here touches a concrete basis or a concrete PP model:
//    alpha_max -- BasisSet::Gaussian::LatticeSum1E::MaxExponent(), documented there as "the GPW
//                 density-grid cutoff floor".
//    alpha_pp  -- BasisSet::SpeciesRadialField_Gaussian::AsGaussians(Z, Short), whose terms carry
//                 alpha = 1/(2 r_loc^2).  A model with no closed-Gaussian short part does not implement
//                 the face; that leaves alpha_pp at 0, which the selector reads as "not measurable" --
//                 NOT as "smooth".
//---------------------------------------------------------------------------------------------------
static qcMesh::XCMeshSharpness GatherSharpness(const Lattice_3D& lat, const BasisSet::Real_BS& mol,
                                               const Pseudopotential::PseudopotentialSet& pp, bool imposed)
{
    qcMesh::XCMeshSharpness s;
    s.cellEdge = lat.GetUnitCell().GetMaximumCellEdge();
    s.nAtoms   = int(lat.GetUnitCell().GetNumAtoms());
    s.imposed  = imposed;   // the RESOLVED value: CP2K_COMPAT can veto it, and the grid cost follows
    // The basis's own answer (V-BS-RASTER): the sharpest function exponent read off its density's resolution, as the ladder reads it --
    // no cast to a Gaussian face, and every block counted (it was the first block that answered).
    {
        const double eps=BasisSet::Gaussian::GpwRasterLadder::kResolutionEps;
        s.alphaMax = BasisSet::Gaussian::ExponentAtCutoff(mol.MakeResolution(eps).cutoff, eps);
    }
    for (const auto& [element, valence] : pp.Species(*lat.GetStructure()))   // the structure's elements, never the whole set
    {
        const int Z = int(thePeriodicTable().GetZ(element));
        const Pseudopotential::HGH_LocalPotential loc = Pseudopotential::GetGTH(element,pp.Functional(),valence).local;
        const auto* g = static_cast<const BasisSet::SpeciesRadialField_Gaussian*>(&loc);
        for (const auto& t : g->AsGaussians(Z, BasisSet::FieldRange::Short)) s.alphaPP = std::max(s.alphaPP, t.alpha);
    }
    return s;
}

//---------------------------------------------------------------------------------------------------
//  THE MAGNETIC ORDERING IS STATED (OpenWork §3 "The magnetic decoration is DERIVED FROM THE SEED", 2026-10-26).
//  The message names the sites, not a derived guess: the ordering is the caller's chemistry.
void RequireStatedMagneticOrdering(const SolidCalcOptions& o, const Structure& st)
{
    std::string sites;
    st.ForEachSite([&](int Z, const rvec3_t&, bool){ sites+=(sites.empty()?"":" ")+thePeriodicTable().GetSymbol(Z); });
    const size_t nAtoms=st.GetNumAtoms();
    if (!o.siteSpins.empty())
    {
        if (o.siteSpins.size()!=nAtoms)
            throw std::invalid_argument("SolidCalcOptions: siteSpins has "+std::to_string(o.siteSpins.size())
                                        +" entries for "+std::to_string(nAtoms)+" sites ("+sites+")");
        for (int s : o.siteSpins)
            if (s<-1 || s>1) throw std::invalid_argument("SolidCalcOptions: siteSpins entries are -1, 0 or +1, got "+std::to_string(s));
        return;
    }
    // R2.5 (doc/SCFRobustnessPlan.md P2, 2026-10-30): an IMPOSED polarized run states it too, +U or not -- the imposition IS
    // the decoration (the Shubnikov group), and deriving it from the seed let SAD impose a ferrimagnet on MnO AFM-II.
    const bool polarized = o.multiplicity>=1;
    if (polarized && (!o.hubbard.empty() || o.imposeSymmetry) && !o.greyImposition)
        throw std::invalid_argument("SolidCalcOptions: a polarized run with a Hubbard manifold or an imposed symmetry must STATE its"
                                    " magnetic ordering -- siteSpins, one of +1/-1/0 per site ("+sites+"); all zeros for a"
                                    " non-magnetic cell.  It is no longer derived from the seed (doc/OpenWork.md §3, R2.5)");
}

//  THE CHARGE ORDER IS STATED OR DERIVED (OpenWork §3 "Spinel charge order is CHOSEN by atom index", 2026-10-30).
void RequireConsistentSiteCharges(const SolidCalcOptions& o, const Structure& st)
{
    if (o.siteCharges.empty()) return;
    std::string sites;
    st.ForEachSite([&](int Z, const rvec3_t&, bool){ sites+=(sites.empty()?"":" ")+thePeriodicTable().GetSymbol(Z); });
    if (o.siteCharges.size()!=st.GetNumAtoms())
        throw std::invalid_argument("SolidCalcOptions: siteCharges has "+std::to_string(o.siteCharges.size())
                                    +" entries for "+std::to_string(st.GetNumAtoms())+" sites ("+sites+")");
    int sum=0; for (int q : o.siteCharges) sum+=q;
    if (sum!=0) throw std::invalid_argument("SolidCalcOptions: siteCharges sum to "+std::to_string(sum)+", not 0 -- formal charges conserve the"
                                            " cell's charge ("+sites+")");
    if (o.seed!=ChargeDensity::SeedStrategy::IonicSAD && o.seed!=ChargeDensity::SeedStrategy::Default)   // Default: the Auto rule takes IonicSAD
        throw std::invalid_argument("SolidCalcOptions: siteCharges is stated but the seed is not IonicSAD -- only the ionic seed reads a"
                                    " charge order, so nothing would consume it");
}

//  THE SPIN BOOKKEEPING (decision 2).  0 = minimal spin (ζ = 0 for an even count, a doublet for an odd one); > 1 = 2S+1 as stated;
//  1 = the EXPLICIT two-channel singlet, nUp = nDown, which an odd count cannot hold.  Until 2026-11-01 1 mapped to nel%2 like 0, so
//  the parity check could never fire on it: bcc Li (1 e/cell) ran a fully spin-up doublet, 53 mHa off, and reported CONVERGED
//  (OpenWork §3).  Under spinsShareFermi the counts are only the split the occupations start from -- one Fermi level decides them --
//  so an odd count keeps the nel%2 split there.
//---------------------------------------------------------------------------------------------------
//  ROUTES THE CODE DOES NOT REALIZE YET, refused from the options alone -- so BEFORE the basis is built (priority 3: refuse
//  before spending).  The words are the deck's (xc.fit, basis.ecut), with the C++ field beside them.
static void RefuseUnrealizedRoutes(const SolidCalcOptions& o)
{
    if (o.vxcFit==Hamiltonian::VxcFit::Auxiliary)
        policy::Refuse("VxcProductSpace", "xc.fit Auxiliary (the Dunlap auxiliary fit) on a periodic cell",
                       "the auxiliary Vxc fit is a molecular / atom fit with no periodic realization yet, an open cell of the Vxc "
                       "product space (doc/Pins.md pin 30) -- state xc.fit Delta or PlaneWave");
    if (o.orbitalBasis==SolidCalcOptions::OrbitalBasis::PlaneWave && !(o.planeWaveEcut>0.0))
        policy::Refuse("PlaneWaveEcut", "basis.family planeWave (SolidCalcOptions::orbitalBasis) with no cutoff",
                       "a plane-wave basis is its cutoff -- state basis.ecut (planeWaveEcut) > 0 Ha");
}

std::pair<int,int> SpinChannelCounts(const SolidCalcOptions& o, int nel)
{
    const bool singlet = o.multiplicity==1 && !o.spinsShareFermi;
    if (singlet && nel%2!=0)
        throw std::invalid_argument("SolidCalcOptions: multiplicity 1 is the two-channel SINGLET (nUp = nDown), impossible with an odd Nelec "
                                    +std::to_string(nel)+" -- state multiplicity 2 (a doublet), 0 (minimal spin), or spinsShareFermi true so"
                                    " one Fermi level sets the channels (a metal)");
    const int twoS = o.multiplicity>1 ? o.multiplicity-1 : (singlet ? 0 : nel%2);
    if ((nel-twoS)%2!=0 || twoS>nel)
        throw std::invalid_argument("SolidCalcOptions: multiplicity "+std::to_string(o.multiplicity)+" parity disagrees with Nelec "
                                    +std::to_string(nel));
    return {(nel+twoS)/2, (nel-twoS)/2};
}

//---------------------------------------------------------------------------------------------------
//  THE RESOLVE PHASE.  What used to be the decision half of SolidCalculation::BuildBasis, in the order it decided them; the
//  build half (the Bloch basis, its report, the electron configuration) is SolidCalculation::BuildPoint.
namespace
{
//! THE STRAIN-ALLOWANCE POLICY (the cell half's frozen rule): the reference lattice with its sizing box set to (1+s) A_ref, reported
//! under policy.  s = 0 returns the lattice untouched (every size decided on the cell itself, bit for bit as before).
Lattice_3D SizedForStrain(const Lattice_3D& lat, const SolidCalcOptions& opts)
{
    const double s=opts.cellStrainAllowance;
    if (!(s>=0.0)) throw std::invalid_argument("SolidCalcOptions: cellStrainAllowance must be >= 0 (a linear fraction, 0.03 = 3%)");
    if (s==0.0) return lat;
    if (s>0.25) throw std::invalid_argument("SolidCalcOptions: cellStrainAllowance "+std::to_string(s)+" > 0.25: the sets sized for a cell that "
                                            "much larger cost (1+s)^3 in every grid -- relax in stages, re-resolving the setup at the new cell");
    UnitCell cell=lat.GetUnitCell();
    cell.SetSizingCell((1.0+s)*cell.GetCellMatrix());
    std::ostringstream why;
    why<<"every G set, raster and mesh count is decided on the reference cell grown by "<<100.0*s<<"% and held as the cell strains (the frozen "
         "rule; CP2K's CELL_REF): E(cell) stays smooth for any cell whose largest principal stretch is <= "<<1.0+s<<", and the reference point "
         "runs on sets ~"<<std::setprecision(3)<<std::pow(1.0+s, 3.0)<<"x the fixed-cell size";
    policy::Report("SolidSetup.CellStrainAllowance", "sizing box (1+"+std::to_string(s)+") A_ref", policy::Verdict::Choose, why.str());
    return Lattice_3D(cell, lat.GetLimits());
}
}

std::shared_ptr<const SolidSetup> SolidSetup::Resolve(const Lattice_3D& lat0, std::shared_ptr<const BasisSet::Real_BS> mol,
                                                      const SolidCalcOptions& opts,
                                                      const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    qchem::SetRunPolicy(opts.policy);   // D-ENV 6b: the policy is part of the run's options, installed ONCE here, before any factory consults it
    RefuseUnrealizedRoutes(opts);
    return Decide(SizedForStrain(lat0, opts), std::move(mol), opts, acc);
}

std::shared_ptr<const SolidSetup> SolidSetup::Resolve(const Lattice_3D& lat0, const GaussianBasisRecipe& recipe,
                                                      const SolidCalcOptions& opts,
                                                      const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    qchem::SetRunPolicy(opts.policy);   // before the vet loop builds its Bloch blocks, as before every other factory
    RefuseUnrealizedRoutes(opts);      // before the basis is built: nothing is spent on a run that cannot go
    const Lattice_3D lat=SizedForStrain(lat0, opts);   // the vet loop's Bloch blocks are sized as every later point's will be
    const auto pp=Pseudopotential::PseudopotentialSet::Get(opts.pp).With(opts.species);
    // THE BASIS, made over a structure the setup will hold (Decide's setup keeps its own copy of the same cell, and the basis is
    // built over THAT copy below so it cannot outlive its atoms).
    GaussianBasisRecipe rec=recipe;
    rec.tolerances=opts.tolerances;   // the options are the one source of the run's tolerances; the recipe binds them on every build
    rec.tolerances.relCutoffSafety=opts.accuracy.relCutoffSafety;   // the pair rule's stiffness is a calibration: Decide writes the same value
    GaussianBasisRecipe frozen=rec;
    const std::shared_ptr<const Structure> cell=lat.GetStructure();
    std::shared_ptr<const BasisSet::Real_BS> mol;
    if (recipe.autoTrim)
    {
        if (!recipe.trim.empty())
            throw std::invalid_argument("GaussianBasisRecipe: autoTrim and a stated trim are exclusive -- state one");
        if (opts.orthoTol<=0.0)
            throw std::invalid_argument("GaussianBasisRecipe: autoTrim needs an orthoTol > 0 to trim at -- at 0 it would drop nothing");
        // THE VET-STAGE TRIM, ONCE, at the reference geometry (pin 22); frozen below so a geometry step rebuilds the SAME shells.
        auto make=[&](const BasisSet::Gaussian::ShellTrim& t){ GaussianBasisRecipe r=rec; r.trim=t; return r.Make(*cell, pp); };
        BasisSet::Lattice::VetTrimResult r;
        {
            qchem::report::Scope here("basisTrim");   // the vet stage files its passes and its verdict one level in
            r=BasisSet::Lattice::VetStageTrim(lat, make, {.images=opts.images, .kShift=opts.kShift, .tol=rec.tolerances}, opts.orthoTol);   // the tolerances make() binds
        }
        frozen.trim=std::move(r.trim);
        frozen.autoTrim=false;
        mol=std::move(r.mol);
    }
    else
    {
        if (!recipe.trim.empty()) { std::cout<<"[basis trim] STATED (no vet loop): "; recipe.trim.Write(std::cout); std::cout<<std::endl; }
        mol=rec.Make(*cell, pp);
    }
    std::shared_ptr<SolidSetup> s=Decide(lat, std::move(mol), opts, acc);
    s->cell  =cell;     // the structure the basis was built over: the setup holds it as long as it holds the basis
    s->recipe=std::move(frozen);
    return s;
}

std::shared_ptr<SolidSetup> SolidSetup::Decide(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                               const SolidCalcOptions& opts,
                                               const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    std::shared_ptr<SolidSetup> s(new SolidSetup(Pseudopotential::PseudopotentialSet::Get(opts.pp).With(opts.species), lat));
    s->opts=opts;
    s->acc =acc;
    s->mol =std::move(mol);
    const std::shared_ptr<const Structure> st=s->cell;
    RequireStatedMagneticOrdering(opts, *st);
    RequireConsistentSiteCharges(opts, *st);
    // THE ELECTRON COUNT IS THE SPECIES' (2026-10-31): sum of the run's per-atom valence minus the STATED net charge.  Derived when
    // unstated; a stated count that disagrees is REFUSED, naming both -- the class of the silent +8 (a stated Na q9 run with the
    // material's q1 count), where nothing downstream could tell a mis-derived count from an intended charge.
    {
        int sumQ=0;
        lat.GetUnitCell().ForEachSite([&](int Z, const rvec3_t&, bool){ sumQ+=s->pp.Valence(thePeriodicTable().GetSymbol(Z)); });
        const int derived=sumQ-opts.charge;
        if (opts.Nelec==0) s->opts.Nelec=derived;
        else if (opts.Nelec!=derived)
            throw std::invalid_argument("SolidCalcOptions: Nelec "+std::to_string(opts.Nelec)+" disagrees with the species: the cell's valence "
                                        "sums to "+std::to_string(sumQ)+(opts.charge!=0 ? " and the stated charge is "+std::to_string(opts.charge) : "")
                                        +", so a neutral"+(opts.charge!=0 ? "-at-that-charge" : "")+" cell holds "+std::to_string(derived)
                                        +".  A charged cell is legal only when STATED (SolidCalcOptions::charge); leave Nelec 0 to derive it");
    }
    s->seedKey=SeedKeyOf(s->pp, *st);
    s->seedKey.siteCharges=opts.siteCharges;
    // THE MAGNETIC ORDERING HAS ONE SOURCE (controller ruling 2026-11-12, OpenWork §3 "One magnetic ordering has TWO sources of truth"):
    // a STATED siteSpins is the run's ordering, and the seed reads it as the imposition does; the structure's per-atom spins are only the
    // default when siteSpins is absent.  Before, an inline cell without atom spins seeded every site +m, which the imposed AFM Shubnikov
    // group projected to m = 0 -- silently.  A stated ordering whose flips differ from the structure's is reported as the override it is.
    s->seedKey.siteSpins=opts.siteSpins;
    if (!opts.siteSpins.empty())
    {
        std::string fromStructure, differs;
        size_t i=0;
        st->ForEachSite([&](int Z, const rvec3_t&, bool spinFlip)
        {
            const int stated=opts.siteSpins[i++];
            fromStructure+=spinFlip ? " -" : " +";
            if ((stated<0)!=spinFlip) differs+=" "+thePeriodicTable().GetSymbol(Z)+std::to_string(i)+(stated<0 ? "(-)" : stated>0 ? "(+)" : "(0)");
        });
        if (!differs.empty())
        {
            std::string stated;
            for (int v : opts.siteSpins) stated+=v<0 ? " -" : v>0 ? " +" : " 0";
            const std::string decided="seed spins: the stated siteSpins ["+stated+" ] override the structure's flips ["+fromStructure+" ] at"+differs;
            const std::string reason="a stated magnetic ordering is the run's ONE source: the seed reads it, as the imposition does";
            policy::Report("SeedSpins", decided, policy::Verdict::Choose, reason);
            std::cout << "[policy] " << decided << " -- " << reason << std::endl;
        }
    }
    // THE SEED: a stated one is honoured; Default takes THE Auto rule both setups share (qchem.SeedPolicy) -- a lattice seeds from the
    // valence library's superposition at its formal charges (stated, else derived: its Madelung field holds them), and a species the
    // library lacks takes the core guess with a WARN (it threw at seed time until 2026-11-12).
    s->opts.seed=ResolveSeed(opts.seed, *st, {.densityFock = true, .periodic = st->HasTranslationSymmetry(), .chargesStated = !opts.siteCharges.empty(),
        .covered = [key=s->seedKey](int Z){ return ChargeDensity::HasAtomicDensity(Z, key.functional, "atomic_valence_densities.json", -1, key.Q(Z)); }});

    // THE RESOLVED IMPOSITION (doc/OpenWork.md N5, user 2026-08-26).  The caller's flag AND the policy's
    // permission: CP2K parity forbids the capability outright, because CP2K does no symmetry work at all
    // (see RunPolicy::SymmetryImposition for the evidence), and every banked recipe asks for it -- so the
    // veto has to live here rather than in each recipe.  ONE name from here on; nothing below reads
    // opts.imposeSymmetry again, or the two would drift.
    s->imposed = opts.imposeSymmetry && theRunPolicy().SymmetryImposition();
    // THE MAGNETIC DECORATION (S3).  An imposition is only SHUBNIKOV if the factory is told which sites
    // carry which spin; without it an "imposed AFM" run star-averages under the SPATIAL group and the
    // order is erased.  ALWAYS STATED (siteSpins) on an imposed polarized run -- RequireStatedMagneticOrdering above
    // refused the run otherwise (R2.5, 2026-10-30: the seed-derived decoration is gone; with SeedStrategy::SAD it
    // labelled neutral O (2p4) magnetic and MnO AFM-II imposed a ferrimagnet, 12 Shubnikov ops, not 24).
    //   greyImposition -> {} DELIBERATELY: that arm is the erasure NEGATIVE CONTROL
    if (s->imposed && !opts.greyImposition) s->siteSpins = opts.siteSpins;   // {} on an unpolarized run that states none
    s->decoration=s->siteSpins;
    if (s->decoration.empty()) s->decoration.assign(lat.GetUnitCell().GetNumAtoms(), 0);

    s->planeWave = opts.orbitalBasis==SolidCalcOptions::OrbitalBasis::PlaneWave;
    // THE RASTER DEMANDS (C2 step 3, doc/GridEpsilonPlan.md §1): eps = 0 states them (the tolerances as given); a stated calibration
    // derives them here; otherwise the run's own pilot audit sizes them after its SCF (rastersResolved false: the pilot starts on the
    // ball's own rasters).  Stating a demand beside a target would be two answers to one question.
    {
        const auto& a=opts.accuracy;
        const auto& t=opts.tolerances;
        if (!(a.eps>=0.0)) throw std::invalid_argument("SolidCalcOptions: accuracy eps must be >= 0 (0 = no target: the stated rasters)");
        const bool stated = t.rasterX!=0.0 || t.vlocRasterX!=0.0 || t.xcRasterX!=0.0;
        if (stated && a.eps>0.0)
            throw std::invalid_argument("SolidCalcOptions: tolerances state a raster demand (rasterX / vlocRasterX / xcRasterX) beside an accuracy "
                                        "target -- the target derives them; state accuracy 0 to give the rasters yourself");
        if (!s->planeWave && a.eps>0.0)
        {
            if (a.eta) {s->opts.tolerances=a.Stated(t); s->rasterSource=SolidSetup::RasterSource::Calibration;}
            else       {s->rastersResolved=false;       s->rasterSource=SolidSetup::RasterSource::Pilot;}
            // V_loc long's routing tolerance from the same target (a STATED vlocEps is honoured, like a stated densityEcut)
            if (s->opts.tolerances.vlocEps==0.0)
            {
                if (!(a.vlocRoutingEta>0.0)) throw std::invalid_argument("SolidCalcOptions: accuracy vlocRoutingEta (delivered / bound) must be > 0");
                s->opts.tolerances.vlocEps=a.VlocEps();
                s->opts.tolerances.vlocEpsRigorous=a.RigorousVlocEps();   // the two pieces, report-only
                s->opts.tolerances.vlocRoutingEta=a.vlocRoutingEta;
            }
        }
        if (!(a.rungGateC>=0.0)) throw std::invalid_argument("SolidCalcOptions: accuracy rungGateC must be >= 0 (0 = no completion rung)");
        s->opts.tolerances.rungGateC=a.rungGateC;   // with or without a target: the ladder's gate is not derived from eps
        if (!(a.relCutoffSafety>0.0)) throw std::invalid_argument("SolidCalcOptions: accuracy relCutoffSafety (the pair->level rule's stiffness) must be > 0");
        s->opts.tolerances.relCutoffSafety=a.relCutoffSafety;   // nor is the pair rule's stiffness
    }
    if (s->planeWave)
    {
        // THE PLANE-WAVE FAMILY (OpenWork §2 PW-1).  The factory builds ONE Γ block and knows no symmetry, so
        // anything else is refused here rather than silently run as Γ / unimposed.
        assert(opts.planeWaveEcut>0.0 && "RefuseUnrealizedRoutes refused a plane-wave run with no cutoff");
        if (lat.MakeKMesh(opts.kShift).size()!=1 || opts.kShift!=rvec3_t(0,0,0))
            policy::Refuse("PlaneWaveGammaOnly", "basis.family planeWave off Γ (a k-mesh or a kShift)",
                           "the plane-wave factory builds one Γ block so far -- state kmesh [1,1,1] and no kShift");
        // A MAGNETIC imposition is the (D_up, D_dn) pair projection on the one Γ block (PW_BasisSet, OpenWork §2 PW-5).  A
        // grey / spatial imposition would need the per-channel density star-average the plane-wave route does not have.
        const bool magnetic=std::any_of(s->siteSpins.begin(), s->siteSpins.end(), [](int v){return v!=0;});
        if (s->imposed && !magnetic)
            policy::Refuse("PlaneWaveMagneticOnly", "basis.family planeWave with an imposed grey / spatial group",
                           "the plane-wave route imposes only a MAGNETIC (Shubnikov) group so far -- state solid.siteSpins, or "
                           "solid.imposeSymmetry false");
        // Plane waves ARE orthonormal: the default Cholesky would factor an identity S into two dense identity
        // factors per spin block.  A stated non-default ortho is obeyed (an A/B door).
        if (s->opts.ortho==qchem::Cholesky) s->opts.ortho=qchem::Orthonormal;
    }

    // DECISION 1 -- the XC quadrature.  Resolve Auto HERE, once, from facts about the run.  Downstream
    // consumers compare ==Becke, so an unresolved Auto would silently read as Uniform; resolving it at the
    // point the spec enters the Hamiltonian is what makes that impossible.
    // The XC-mesh route is a DECLARED CP2K deviation (RunPolicy::BeckeXC), resolved here beside the
    // imposition for the same reason: the policy belongs to the run, the sizing belongs to qcMesh.
    // A plane-wave run's Auto is its own FFT raster (QE's route): the sharpness gather prices Gaussians, and a
    // delta fit would collocate every PW at every mesh point.  A STATED grid is obeyed, as for GPW.
    // GEOMETRY-INDEPENDENT by construction: the sharpness is exponents (basis, PP) and the cell, never positions.
    if (s->planeWave && opts.xcMesh.cellKind==qcMesh::UnitCellKind::Auto)
        s->xcMesh = qcMesh::MeshParams{.cellKind=qcMesh::UnitCellKind::Uniform};
    else
        s->xcMesh = qcMesh::ResolveXCMesh(opts.xcMesh, GatherSharpness(lat, *s->mol, s->pp, s->imposed),
                                          theRunPolicy().BeckeXC());
    // The Vxc product space's admission (pin 30), asked of the Hamiltonian's rule HERE, before the Hamiltonian is built: the
    // Factory would throw the same reason later, after the basis and its fits.
    if (const auto a=Hamiltonian::AdmitVxc(*st, opts.vxcFit, s->xcMesh); !a)
        policy::Refuse("AdmitVxc", "the stated Vxc fit and XC mesh on this cell", a.Error());
    {   // the rho axis, on the fit the Hamiltonian RESOLVES (Ham_PW_DFT: Auto -> Delta on a Becke mesh, else the plane-wave fit)
        using Hamiltonian::VxcFit;
        const bool becke=s->xcMesh.cellKind==qcMesh::UnitCellKind::Becke;
        const VxcFit fit= opts.vxcFit==VxcFit::Auto ? (becke ? VxcFit::Delta : VxcFit::PlaneWave) : opts.vxcFit;
        ReportSolidXCRho(fit, opts.xcRho, s->planeWave);
    }

    // DECISION 2 -- the spin bookkeeping.  multiplicity -> (nUp,nDown), with the parity check that catches
    // a singlet asked of an odd electron count BEFORE integer division silently empties a channel.
    const int nel  = s->opts.Nelec;   // the RESOLVED count (above)
    s->polarized = opts.multiplicity>=1;
    std::tie(s->nUp, s->nDown)=SpinChannelCounts(opts, nel);
    // A partial eigensolve must leave EMPTY bands above the fuller channel (smearing and the gap need them).
    if (opts.numBands>0 && opts.numBands<=size_t(s->nUp))
        throw std::invalid_argument("SolidCalcOptions: numBands "+std::to_string(opts.numBands)+" must exceed the larger spin channel's "
                                    +std::to_string(s->nUp)+" electrons");
    return s;
}

std::shared_ptr<const BasisSet::Real_BS> GaussianBasisRecipe::Make(const Structure& cell,
                                                                  const Pseudopotential::PseudopotentialSet& pp) const
{
    using namespace BasisSet::Gaussian;
    // D-BASIS-PP: the basis file's declared pseudopotential q per element must be the q this run's terms use.
    const auto species=pp.Species(cell);
    if (file.empty()) CheckPPProvenance(data, species); else CheckPPProvenance(file, species);
    {   // ... and the file must hold every element (the policy row names the bases that do), before anything is read from it
        const std::string name=file.empty() ? BasisName(data) : file;   // the deck's own word: the registered name, or the stated file
        RequireBasisCovers(file.empty() ? name+".bsd" : file, name, cell);
    }
    std::unique_ptr<BasisSet::Real_BS> built(file.empty() ? Factory(data, &cell, Engine::MnD, Angular::Cartesian, trim)
                                                          : FactoryFromFile(file, &cell, Engine::MnD, Angular::Cartesian, trim));
    // D-ENV-TOL option A: bind the run's pair-loop tolerances while the basis is still OURS and mutable -- before a view wraps it and
    // before anything computes from it (a MnD Cartesian basis can take them; anything else is a build defect, not an input).
    auto* bind=dynamic_cast<ToleranceBinding*>(built.get());
    if (!bind) throw std::logic_error("GaussianBasisRecipe::Make: the built basis cannot bind GPW tolerances (not a MnD Cartesian set)");
    bind->BindTolerances(tolerances);
    std::shared_ptr<const BasisSet::Real_BS> b(std::move(built));
    if (angular==Angular::Spherical) b=PG_Spherical::MakeSphericalLatticeView(std::move(b));
    return b;
}

std::vector<Symmetry::SymOp> SolidSetup::KeptOps(const Lattice_3D& at) const
{
    if (!imposed) return {Symmetry::SymOp{Matrix3D<double>(1,0,0, 0,1,0, 0,0,1), rvec3_t(0,0,0)}};
    if (!siteSpins.empty()) return at.ShubnikovOps(siteSpins);
    std::vector<Symmetry::SymOp> ops;
    for (const auto& op : at.GetSpaceGroup().Ops()) ops.push_back({op.W, op.tau});
    return ops;
}

//! The recipe with its trim cleared: the basis every what-if question starts from.
static GaussianBasisRecipe Untrimmed(const std::optional<GaussianBasisRecipe>& recipe, bool planeWave, const char* who)
{
    if (planeWave) throw std::invalid_argument(std::string(who)+": a plane-wave setup has no Gaussian shells to trim");
    if (!recipe)   throw std::invalid_argument(std::string(who)+": a setup from a pre-built basis has no recipe to clear the trim of");
    GaussianBasisRecipe bare=*recipe;
    bare.trim={};
    bare.autoTrim=false;
    return bare;
}

std::unique_ptr<BasisSet::Lattice::ShellConditioning> SolidSetup::MakeShellConditioning() const
{
    const GaussianBasisRecipe bare=Untrimmed(recipe, planeWave, "SolidSetup::MakeShellConditioning");
    return std::make_unique<BasisSet::Lattice::ShellConditioning>(reference, bare.Make(*cell, pp),
                                                                  BasisSet::Lattice::GPWParams{.images=opts.images, .kShift=opts.kShift, .tol=bare.tolerances});
}

BasisSet::Gaussian::ShellTrim SolidSetup::MakeVetTrim(double orthoTol) const
{
    const GaussianBasisRecipe bare=Untrimmed(recipe, planeWave, "SolidSetup::MakeVetTrim");
    if (!(orthoTol>0.0)) throw std::invalid_argument("SolidSetup::MakeVetTrim: an orthoTol > 0 to trim at -- at 0 it would drop nothing");
    auto make=[&](const BasisSet::Gaussian::ShellTrim& t){ GaussianBasisRecipe r=bare; r.trim=t; return r.Make(*cell, pp); };
    return BasisSet::Lattice::VetStageTrim(reference, make, {.images=opts.images, .kShift=opts.kShift, .tol=bare.tolerances}, orthoTol).trim;
}

} //namespace qchem

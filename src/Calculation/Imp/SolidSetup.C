// File: Calculation/Imp/SolidSetup.C  The periodic run setup: every decision made ONCE, before a geometry point is built.
module;
#include <algorithm>   // std::max (alpha_pp)
#include <iostream>    // the stated-trim line
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.SolidSetup;

import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;   // GaussianSharpness::MaxExponent -- alpha_max (the NARROW face)
import qchem.BasisSet.Orbital_1E_IBS;          // Real_OIBS (the per-irrep bases to iterate)
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH local PP -> alpha_pp
import qchem.PeriodicTable;                    // thePeriodicTable().GetZ (element symbol -> Z)
import qchem.Mesh.XCPolicy;                    // XCMeshSharpness / ResolveXCMesh (the grid decision)
import qchem.BasisSet.Gaussian.Lattice.SphericalLatticeView;   // MakeSphericalLatticeView (GaussianBasisRecipe::spherical)

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
    for (auto ibs : const_cast<BasisSet::Real_BS&>(mol).Iterate<BasisSet::Real_OIBS>())
        // THE NARROW FACE, since the 2026-09-08 ISP split: this asks for ONE number, so it names the
        // three-method GaussianSharpness capability and not the seventeen-method periodic aggregate it
        // used to cross-cast to.  Nothing here wants a lattice sum, a collocation or a stream fold.
        if (const auto* sh=dynamic_cast<const BasisSet::Gaussian::GaussianSharpness*>(ibs))
            { s.alphaMax = sh->MaxExponent(); break; }
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
    const bool polarized = o.multiplicity>=1;
    if (polarized && !o.hubbard.empty() && !o.greyImposition)
        throw std::invalid_argument("SolidCalcOptions: a polarized run with a Hubbard manifold must STATE its magnetic ordering --"
                                    " siteSpins, one of +1/-1/0 per site ("+sites+"); all zeros for a non-magnetic cell."
                                    "  It is no longer derived from the seed (doc/OpenWork.md §3)");
}

//---------------------------------------------------------------------------------------------------
//  THE RESOLVE PHASE.  What used to be the decision half of SolidCalculation::BuildBasis, in the order it decided them; the
//  build half (the Bloch basis, its report, the electron configuration) is SolidCalculation::BuildPoint.
std::shared_ptr<const SolidSetup> SolidSetup::Resolve(const Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                                      const SolidCalcOptions& opts,
                                                      const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    qchem::SetRunPolicy(opts.policy);   // D-ENV 6b: the policy is part of the run's options, installed ONCE here, before any factory consults it
    return Decide(lat, std::move(mol), opts, acc);
}

std::shared_ptr<const SolidSetup> SolidSetup::Resolve(const Lattice_3D& lat, const GaussianBasisRecipe& recipe,
                                                      const SolidCalcOptions& opts,
                                                      const SCFAccelerators::SolidAcceleratorOptions& acc)
{
    qchem::SetRunPolicy(opts.policy);   // before the vet loop builds its Bloch blocks, as before every other factory
    const auto pp=Pseudopotential::PseudopotentialSet::Get(opts.pp).With(opts.species);
    // THE BASIS, made over a structure the setup will hold (Decide's setup keeps its own copy of the same cell, and the basis is
    // built over THAT copy below so it cannot outlive its atoms).
    GaussianBasisRecipe frozen=recipe;
    const std::shared_ptr<const Structure> cell=lat.GetStructure();
    std::shared_ptr<const BasisSet::Real_BS> mol;
    if (recipe.vet)
    {
        if (!recipe.trim.empty())
            throw std::invalid_argument("GaussianBasisRecipe: vet and a stated trim are exclusive -- state one");
        // THE VET-STAGE TRIM, ONCE, at the reference geometry (pin 22); frozen below so a geometry step rebuilds the SAME shells.
        auto make=[&](const BasisSet::Gaussian::ShellTrim& t){ GaussianBasisRecipe r=recipe; r.trim=t; return r.Make(*cell, pp); };
        auto r=BasisSet::Lattice::VetStageTrim(lat, make, {.images=opts.images, .kShift=opts.kShift}, opts.orthoTol);
        frozen.trim=std::move(r.trim);
        frozen.vet =false;
        mol=std::move(r.mol);
    }
    else
    {
        if (!recipe.trim.empty()) { std::cout<<"[basis trim] STATED (no vet loop): "; recipe.trim.Write(std::cout); std::cout<<std::endl; }
        mol=recipe.Make(*cell, pp);
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
    s->seedKey=SeedKeyOf(s->pp, *st);

    // THE RESOLVED IMPOSITION (doc/OpenWork.md N5, user 2026-08-26).  The caller's flag AND the policy's
    // permission: CP2K parity forbids the capability outright, because CP2K does no symmetry work at all
    // (see RunPolicy::SymmetryImposition for the evidence), and every banked recipe asks for it -- so the
    // veto has to live here rather than in each recipe.  ONE name from here on; nothing below reads
    // opts.imposeSymmetry again, or the two would drift.
    s->imposed = opts.imposeSymmetry && theRunPolicy().SymmetryImposition();
    // THE MAGNETIC DECORATION (S3).  An imposition is only SHUBNIKOV if the factory is told which sites
    // carry which spin; without it an "imposed AFM" run star-averages under the SPATIAL group and the
    // order is erased.  STATED (siteSpins) wherever a Hubbard manifold rides along -- RequireStatedMagneticOrdering
    // above refused the run otherwise.  The derivation from the seed below survives ONLY for an imposed polarized
    // run WITHOUT +U, and it is the defect OpenWork §3 names: with SeedStrategy::SAD it labels neutral O (2p4)
    // magnetic and MnO AFM-II imposes a ferrimagnet (12 Shubnikov ops, not 24).
    //   unpolarized  -> {} (no channels to decorate)
    //   greyImposition -> {} DELIBERATELY: that arm is the erasure NEGATIVE CONTROL
    if (!s->imposed) { /* nothing to decorate: siteSpins are meaningless without an imposition */ }
    else if (!opts.siteSpins.empty() && !opts.greyImposition)
        s->siteSpins = opts.siteSpins;              // STATED: a specific ordering, not the seed's guess
    else if (opts.multiplicity>=1 && !opts.greyImposition)
    {
        const std::vector<int> targets = (opts.seed==qchem::ChargeDensity::SeedStrategy::IonicSAD)
                                       ? qchem::ChargeDensity::IonicSADTargets(st.get(), s->seedKey)
                                       : std::vector<int>{};
        s->siteSpins = qchem::ChargeDensity::MagneticDecoration(st.get(), s->seedKey, targets);
    }
    s->decoration=s->siteSpins;
    if (s->decoration.empty()) s->decoration.assign(lat.GetUnitCell().GetNumAtoms(), 0);

    s->planeWave = opts.orbitalBasis==SolidCalcOptions::OrbitalBasis::PlaneWave;
    if (s->planeWave)
    {
        // THE PLANE-WAVE FAMILY (OpenWork §2 PW-1).  The factory builds ONE Γ block and knows no symmetry, so
        // anything else is refused here rather than silently run as Γ / unimposed.
        if (!(opts.planeWaveEcut>0.0))
            throw std::invalid_argument("SolidCalcOptions: orbitalBasis PlaneWave needs planeWaveEcut > 0 (Ha)");
        if (lat.MakeKMesh(opts.kShift).size()!=1 || opts.kShift!=rvec3_t(0,0,0))
            throw std::invalid_argument("SolidCalcOptions: the plane-wave family is Γ-only so far (k-mesh 1x1x1, no shift)");
        if (s->imposed)
            throw std::invalid_argument("SolidCalcOptions: the plane-wave family has no symmetry imposition yet");
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

    // DECISION 2 -- the spin bookkeeping.  multiplicity -> (nUp,nDown), with the parity check that catches
    // a singlet asked of an odd electron count BEFORE integer division silently empties a channel.
    const int twoS = opts.multiplicity>1 ? opts.multiplicity-1 : opts.Nelec%2;
    s->polarized = opts.multiplicity>=1;
    if ((opts.Nelec-twoS)%2!=0 || twoS>opts.Nelec)
        throw std::runtime_error("SolidCalcOptions: multiplicity "+std::to_string(opts.multiplicity)
                                 +" parity disagrees with Nelec "+std::to_string(opts.Nelec));
    // A partial eigensolve must leave EMPTY bands above the fuller channel (smearing and the gap need them).
    if (opts.numBands>0 && opts.numBands<=size_t((opts.Nelec+twoS)/2))
        throw std::invalid_argument("SolidCalcOptions: numBands "+std::to_string(opts.numBands)+" must exceed the larger spin channel's "
                                    +std::to_string((opts.Nelec+twoS)/2)+" electrons");
    s->nUp  =(opts.Nelec+twoS)/2;
    s->nDown=(opts.Nelec-twoS)/2;
    return s;
}

std::shared_ptr<const BasisSet::Real_BS> GaussianBasisRecipe::Make(const Structure& cell,
                                                                  const Pseudopotential::PseudopotentialSet& pp) const
{
    using namespace BasisSet::Gaussian;
    // D-BASIS-PP: the basis file's declared pseudopotential q per element must be the q this run's terms use.
    const auto species=pp.Species(cell);
    if (file.empty()) CheckPPProvenance(data, species); else CheckPPProvenance(file, species);
    std::shared_ptr<const BasisSet::Real_BS> b(file.empty() ? Factory(data, &cell, Engine::MnD, Angular::Cartesian, trim)
                                                            : FactoryFromFile(file, &cell, Engine::MnD, Angular::Cartesian, trim));
    if (spherical) b=PG_Spherical::MakeSphericalLatticeView(std::move(b));
    return b;
}

} //namespace qchem

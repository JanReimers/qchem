// File: Calculation/Imp/BasisMenu.C  The orbital-basis menu of a structure, composed from the policy layer's own facts.
module;
#include <set>
#include <string>
#include <vector>
module qchem.BasisMenu;
import qchem.BasisCoverage;                              // BasisCovers, ElementsOf: the reader's own coverage test
import qchem.BasisSet.Gaussian.Point.BasisFiles;         // BasisNames, KindOfBasisFile: the file's stated KIND
import qchem.BasisSet.Gaussian.Point.Factory;            // RegisteredBasis, CartesianOnlyElements
import qchem.Symmetry.Group;                             // HasLattice, kAdmissionTolerance: is the structure periodic (the route rule)
import qchem.PeriodicTable;
import qchem.Symmetry.Molecule.PointGroup;               // DetectPointGroup: the molecule's SALC blocking group

namespace qchem
{
namespace
{
//! The radial atom's families, in the deck's spelling (Deck.C kAtomType; RadialAtomBasisNames' unit test holds the two together).
struct RadialFamily { const char* name; const char* note; };
constexpr RadialFamily theRadialFamilies[] = {
    {"Slater",       ""},
    {"Gaussian",     ""},
    {"BSpline6",     ""},
    {"BSpliner6",    ""},
    {"Gaussian_RKB", "four-component (RKB): relativistic"},
    {"Slater_RKB",   "four-component (RKB): relativistic"},
};
constexpr const char* theRadialBlocks = "O(3) (by l, m analytic)";   // the radial atom's blocks: one per l, the m sum analytic

//! The registered Gaussian basis files (fit and test files are not menu entries), as the BasisCoverage refusal suggests them.
std::vector<std::string> RegisteredGaussianNames()
{
    std::vector<std::string> out;
    for (const std::string& name : BasisSet::Gaussian::BasisNames())
    {
        BasisSet::Gaussian::BasisSetData data{};
        if (BasisSet::Gaussian::RegisteredBasis(name, data)) out.push_back(name);
    }
    return out;
}

std::string Symbols(const std::vector<int>& Zs)
{
    std::string out;
    for (int Z : Zs) out += (out.empty() ? "" : ", ") + thePeriodicTable().GetSymbol(Z);
    return out;
}

//! One Gaussian file on route \a route: admitted when it holds a block for every element (the reader's own test) and, on the
//! crystal route, when the file states it is a pseudopotential-valence basis (its KIND line: a crystal runs on a pseudopotential).
//! The file's Cartesian-only elements are its note (their s span lives in the d contaminants: the spherical view would lose it).
BasisOffer GaussianOffer(const std::string& route, const std::string& name, const std::set<int>& elements, const std::string& blocks)
{
    using BasisSet::Gaussian::BasisKind;
    const BasisSet::Gaussian::BasisFileKind kind = BasisSet::Gaussian::KindOfBasisFile(name+".bsd");
    std::vector<int> missing;
    for (int Z : elements) if (!BasisCovers(name+".bsd", Z)) missing.push_back(Z);
    std::string reason = missing.empty() ? std::string() : name+" has no block for "+Symbols(missing);
    if (reason.empty() && route=="crystal" && kind.kind!=BasisKind::PPValence)
        reason = kind.kind==BasisKind::Unstated
            ? name+" does not state its KIND: a crystal runs on a pseudopotential and needs a pseudopotential-valence basis"
            : name+" is an "+BasisSet::Gaussian::KindName(kind.kind)+" basis: a crystal runs on a pseudopotential and needs a "
              "pseudopotential-valence basis";
    std::string note;
    for (const std::string& el : BasisSet::Gaussian::CartesianOnlyElements(name+".bsd"))
        note += (note.empty() ? "Cartesian d only for " : ", ") + el;
    return BasisOffer{.route=route, .key="basis.data="+name, .name=name, .admitted=reason.empty(), .reason=reason,
                      .blockGroup=blocks, .note=note, .optIn=false,
                      .kind=BasisSet::Gaussian::KindName(kind.kind), .ppSet=kind.ppSet};
}
} // anonymous namespace

std::vector<std::string> RadialAtomBasisNames()
{
    std::vector<std::string> out;
    for (const RadialFamily& f : theRadialFamilies) out.push_back(f.name);
    return out;
}

std::vector<BasisOffer> OrbitalBasisOffers(const Structure& st, double symmetryTol)
{
    std::vector<BasisOffer> out;
    const std::set<int> elements = ElementsOf(st);
    // THE ROUTE (deck::KindOf's rule): periodic -> the crystal; one site -> the radial atom; else the molecule.  Periodicity is
    // asked of the structure's Group: does it have a lattice?
    if (st.MakeSymmetry(Symmetry::kAdmissionTolerance)->HasLattice())
    {
        for (const std::string& name : RegisteredGaussianNames())
            out.push_back(GaussianOffer("crystal", name, elements, "T (by k)"));
        out.push_back(BasisOffer{.route="crystal", .key="basis.family=planeWave", .name="plane waves", .admitted=true, .reason="",
                                 .blockGroup="T (by k)", .note="", .optIn=false, .kind="", .ppSet=""});
        return out;
    }
    if (st.GetNumAtoms()==1)
    {
        for (const RadialFamily& f : theRadialFamilies)
            out.push_back(BasisOffer{.route="atom", .key=std::string("atom.type=")+f.name, .name=f.name, .admitted=true, .reason="",
                                     .blockGroup=theRadialBlocks, .note=f.note, .optIn=false, .kind="", .ppSet=""});
        // Q1 (user, 2026-11-21): the one-atom MOLECULE, a route entry offered only when asked; its menu is the molecule's
        out.push_back(BasisOffer{.route="molecule", .key="", .name="3-D Gaussians (one-atom molecule)", .admitted=true, .reason="",
                                 .blockGroup="D2h when imposed, else {e}",
                                 .note="3-D Gaussians: shows broken-symmetry shapes; slower", .optIn=true, .kind="", .ppSet=""});
        return out;
    }
    return MoleculeBasisOffers(st, symmetryTol);
}

std::vector<BasisOffer> MoleculeBasisOffers(const Structure& st, double symmetryTol)
{
    std::vector<BasisOffer> out;
    std::vector<Symmetry::Molecule::SymPoint> pts;
    st.ForEachSite([&](int Z, const rvec3_t& R, bool){ pts.push_back({Z, R}); });
    const auto g = Symmetry::Molecule::DetectPointGroup(pts, symmetryTol);
    const std::string blocks = g.abelian+" when imposed (the abelian subgroup of "+g.symbol+"), else {e}";
    const std::set<int> elements = ElementsOf(st);
    for (const std::string& name : RegisteredGaussianNames())
        out.push_back(GaussianOffer("molecule", name, elements, blocks));
    return out;
}

} // namespace qchem

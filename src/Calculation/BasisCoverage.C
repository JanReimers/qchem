// File: Calculation/BasisCoverage.C  The basis-coverage policy: a basis that lacks an element of the structure is refused at SETUP,
// before anything is built, naming the element, the basis and the registered bases that DO cover it.
//
// The test is the reader's own (Gaussian94Reader::FindAtom), so the pre-flight cannot disagree with the build that would follow --
// it only moves the refusal ahead of the work and says what to do.  A policy of the policy layer (doc/OpenWork.md, the policy row):
// it files a `policy` row, BasisCoverage, before it throws.
module;
#include <algorithm>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.BasisCoverage;
import qchem.BasisSet.Gaussian.Point.BasisFiles;        // BasisFile, BasisNames
import qchem.BasisSet.Gaussian.Point.Factory;           // RegisteredBasis: the orbital bases worth suggesting
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94; // Gaussian94Reader::FindAtom: the build's own test
import qchem.PeriodicTable;
import qchem.Policy;
import qchem.Structure;

export namespace qchem
{

//! Does the basis file \a bsdFile (a bare name in the data directory, or a path) hold a block for element \a Z?
inline bool BasisCovers(const std::string& bsdFile, int Z)
{
    try
    {
        BasisSet::Gaussian::Gaussian94Reader reader(BasisSet::Gaussian::BasisFile(bsdFile));
        reader.FindAtom(Site{Z});
        return true;
    }
    catch (const std::runtime_error&) { return false; }
}

//! The registered orbital bases (fit and test files are not suggested) that hold a block for every element of \a elements, comma
//! separated ("" when none does).
inline std::string BasesCovering(const std::set<int>& elements)
{
    std::string out;
    for (const std::string& name : BasisSet::Gaussian::BasisNames())
    {
        BasisSet::Gaussian::BasisSetData data{};
        if (!BasisSet::Gaussian::RegisteredBasis(name, data)) continue;
        if (std::all_of(elements.begin(), elements.end(), [&](int Z){ return BasisCovers(name+".bsd", Z); }))
            out += (out.empty() ? "" : ", ") + name;
    }
    return out;
}

//! The elements of \a st, by Z.
inline std::set<int> ElementsOf(const Structure& st)
{
    std::set<int> elements;
    st.ForEachSite([&](int Z, const auto&, bool){ elements.insert(Z); });
    return elements;
}

//! Refuse a run that STATES NO BASIS: there is no default (a default that cannot build is a silent trap -- CalcOptions' old "sto-3g"
//! had no file), and a refusal that names the bases covering this structure costs the user one line.  Policy row BasisUnstated.
inline void RequireBasisStated(const std::string& basis, const Structure& st)
{
    if (!basis.empty()) return;
    std::string symbols;
    for (int Z : ElementsOf(st)) symbols += (symbols.empty() ? "" : ", ") + thePeriodicTable().GetSymbol(Z);
    const std::string covering = BasesCovering(ElementsOf(st));
    policy::Refuse("BasisUnstated", "no basis stated",
                   covering.empty() ? "there is no default basis; no registered basis covers "+symbols+" -- state a basis file"
                                    : "there is no default basis; these cover "+symbols+": "+covering);
}

//! Refuse \a st on the basis \a bsdFile (shown as \a basisName, the deck's own word) when some element has no block: a policy row
//! BasisCoverage, then std::invalid_argument naming each missing element and the registered bases that cover all of them.
inline void RequireBasisCovers(const std::string& bsdFile, const std::string& basisName, const Structure& st)
{
    const std::set<int> elements = ElementsOf(st);
    std::vector<int> missing;
    for (int Z : elements) if (!BasisCovers(bsdFile, Z)) missing.push_back(Z);
    if (missing.empty()) return;
    std::string what;
    for (int Z : missing) what += (what.empty() ? "" : ", ") + thePeriodicTable().GetSymbol(Z);
    const std::string alternatives = BasesCovering(elements);
    policy::Refuse("BasisCoverage", "basis "+basisName+" has no block for "+what,
                   alternatives.empty() ? std::string("no registered basis covers every element of this structure -- state a basis "
                                                      "file (basis.file) that does")
                                        : "state a basis that covers every element: "+alternatives);
}

} // namespace qchem

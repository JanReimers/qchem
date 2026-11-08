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

//! Refuse \a st on the basis \a bsdFile (shown as \a basisName, the deck's own word) when some element has no block: a policy row
//! BasisCoverage, then std::invalid_argument naming each missing element and the registered bases that cover all of them.
inline void RequireBasisCovers(const std::string& bsdFile, const std::string& basisName, const Structure& st)
{
    std::set<int> elements;
    st.ForEachSite([&](int Z, const auto&, bool){ elements.insert(Z); });
    std::vector<int> missing;
    for (int Z : elements) if (!BasisCovers(bsdFile, Z)) missing.push_back(Z);
    if (missing.empty()) return;
    std::string what, alternatives;
    for (int Z : missing) what += (what.empty() ? "" : ", ") + thePeriodicTable().GetSymbol(Z);
    for (const std::string& name : BasisSet::Gaussian::BasisNames())
    {
        BasisSet::Gaussian::BasisSetData data{};
        if (!BasisSet::Gaussian::RegisteredBasis(name, data)) continue;   // a fit or test file is no orbital basis to suggest
        if (std::all_of(elements.begin(), elements.end(), [&](int Z){ return BasisCovers(name+".bsd", Z); }))
            alternatives += (alternatives.empty() ? "" : ", ") + name;
    }
    policy::Refuse("BasisCoverage", "basis "+basisName+" has no block for "+what,
                   alternatives.empty() ? std::string("no registered basis covers every element of this structure -- state a basis "
                                                      "file (basis.file) that does")
                                        : "state a basis that covers every element: "+alternatives);
}

} // namespace qchem

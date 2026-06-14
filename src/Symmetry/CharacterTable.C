// File: Symmetry/CharacterTable.C  Character tables of the abelian (D2h-family) point groups.
//
// Stage 3 of the molecular-symmetry plan.  SCF blocking uses only the 8 abelian point groups
// (D2h and its subgroups); their irreps are all 1-dimensional with characters +/-1, so each
// table is a small +/-1 matrix.  The irrep names are the Mulliken labels that will decorate
// the orbital eigenvalue / occupation tables.  Operations are identified by canonical tags
// (E, C2z, C2y, C2x, i, sxy, sxz, syz, sh) so the geometric builder can produce a concrete
// SymOp per column and the SALC projector can pair chi(g) with the operation rep M(g).
module;
#include <string>
#include <vector>
export module qchem.Symmetry.CharacterTable;

export namespace Symmetry
{

struct CharacterTable
{
    std::string                   symbol;   // abelian group, e.g. "C2v"
    std::vector<std::string>      opTags;   // canonical operation tags, one per column
    std::vector<std::string>      irreps;   // Mulliken labels, one per row
    std::vector<std::vector<int>> chi;      // chi[irrep][op] in {+1,-1}; size nIrreps x order

    size_t order()   const { return opTags.size(); }   // group order h (abelian: = #classes)
    size_t nIrreps() const { return irreps.size(); }
};

// The character table of one of the 8 abelian point groups: C1, Ci, Cs, C2, C2h, C2v, D2, D2h.
CharacterTable AbelianCharacterTable(const std::string& symbol);

} //namespace

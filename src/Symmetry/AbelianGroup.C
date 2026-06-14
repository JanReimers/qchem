// File: Symmetry/AbelianGroup.C  The concrete abelian point group used for SCF blocking.
//
// Stage 3a-ii: detect the full group, take its abelian-subgroup character table, and build a
// concrete SymOp for each operation tag in the molecule's own axis frame (principal axis = z,
// plus secondary C2 axes / mirror normals as x,y).  Every constructed operation is verified
// against the geometry.  The result (character table + matching operations) is exactly what
// the SALC projector consumes: chi(g) paired with the operation rep M(g).
module;
#include <vector>
export module qchem.Symmetry.AbelianGroup;
export import qchem.Symmetry.PointGroup;       // SymPoint, SymOp, DetectPointGroup, finders
export import qchem.Symmetry.CharacterTable;   // CharacterTable

export namespace Symmetry
{

struct AbelianGroup
{
    CharacterTable     table;   // symbol, opTags, irreps (Mulliken labels), characters
    std::vector<SymOp> ops;     // concrete operations aligned to the molecule, one per opTag
};

// Build the concrete abelian point group of a molecule (the subgroup used for SCF blocking).
AbelianGroup BuildAbelianGroup(const std::vector<SymPoint>& pts, double tol);

} //namespace

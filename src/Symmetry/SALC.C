// File: Symmetry/SALC.C  Symmetry-adapted linear combinations (the transform O).
//
// Stage 3b.  For each irrep Gamma of the abelian group, the projection operator
//   P^Gamma = (1/h) sum_g chi^Gamma(g) M(g)
// (M(g) = the AO operation rep from BuildOperationRep) projects the AO space onto the irrep-
// Gamma subspace.  An orthonormal basis of P^Gamma's column space gives that irrep's SALC
// columns; concatenating the blocks over all irreps yields the full transform O (nAO x nAO),
// block-structured by irrep.  O block-diagonalizes every symmetry-commuting matrix (S, T, V,
// Fock) into per-irrep blocks; within-block Lowdin orthonormalization stays downstream
// (LASolver), exactly as for atoms.  Each column carries its Mulliken irrep label, which will
// decorate the orbital eigenvalue / occupation tables.
module;
#include <vector>
#include <string>
export module qchem.Symmetry.SALC;
export import qchem.Symmetry.AbelianGroup;   // AbelianGroup (SymOp, CharacterTable)
export import qchem.Symmetry.OperationRep;   // AoShell, BuildOperationRep, rmat_t (re-exports Cart/Sph reps)

export namespace qchem::Symmetry
{

struct SALCs
{
    rmat_t                   O;          // nAO x nAO, columns = SALCs grouped by irrep
    std::vector<std::string> irrep;      // Mulliken irrep label of each column (size nAO)
    std::vector<size_t>      blockStart; // first column of each irrep block (size nIrreps+1)
};

// Build the SALC transform O for an AO basis (described by `shells`) under the abelian group g.
SALCs BuildSALCs(const std::vector<AoShell>& shells, const AbelianGroup& g,
                 const rvec3_t& origin, double tol);

} //namespace

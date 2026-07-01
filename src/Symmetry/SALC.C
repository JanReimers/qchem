//! \file
//! \brief Symmetry-adapted linear combinations (the transform \f$O\f$).
//!
//! For each irrep \f$\Gamma\f$ of the abelian group, the projection operator
//! \f[ P^\Gamma = \frac{1}{h}\sum_g \chi^\Gamma(g)\,M(g) \f]
//! (\f$M(g)\f$ = the AO operation rep from \c BuildOperationRep) projects the AO space onto the
//! \f$\Gamma\f$-subspace.  An orthonormal basis of \f$\operatorname{col}(P^\Gamma)\f$ gives that irrep's
//! SALC columns; concatenating the blocks over all irreps yields the full transform \f$O\f$
//! (\f$n_{AO}\times n_{AO}\f$), block-structured by irrep.  \f$O\f$ block-diagonalizes every
//! symmetry-commuting matrix (\f$S,T,V,\f$ Fock) into per-irrep blocks; within-block Löwdin
//! orthonormalization stays downstream (LASolver), exactly as for atoms.  Each column carries its
//! Mulliken irrep label, which decorates the orbital eigenvalue / occupation tables.
module;
#include <vector>
#include <string>
export module qchem.Symmetry.SALC;
export import qchem.Symmetry.AbelianGroup;   // AbelianGroup (SymOp, CharacterTable)
export import qchem.Symmetry.OperationRep;   // AoShell, BuildOperationRep, rmat_t (depends only on ShellRep)

export namespace qchem::Symmetry
{

//! The SALC transform and its irrep bookkeeping.
struct SALCs
{
    rmat_t                   O;          //!< \f$n_{AO}\times n_{AO}\f$; columns = SALCs grouped by irrep
    std::vector<std::string> irrep;      //!< Mulliken irrep label of each column (size \f$n_{AO}\f$)
    std::vector<size_t>      blockStart; //!< first column of each irrep block (size \f$n_{irreps}+1\f$)
};

//! \brief Build the SALC transform \f$O\f$ for an AO basis (described by \a shells) under the abelian
//! group \a g, about \a origin (the point set centroid), matching centres to a tolerance \a tol.
SALCs BuildSALCs(const std::vector<AoShell>& shells, const AbelianGroup& g,
                 const rvec3_t& origin, double tol);

} //namespace

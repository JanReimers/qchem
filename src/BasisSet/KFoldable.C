// File: BasisSet/KFoldable.C  An orbital block that can share its real-space work with the other k-blocks of its lattice basis (the k-fold).
module;
#include <array>
#include <cstdint>
#include <vector>
export module qchem.BasisSet.KFoldable;
export import qchem.Types;   // ivec3_t, dcmplx, chmat_t

export namespace qchem::BasisSet
{

//! \brief A k-BLOCK OF ONE LATTICE BASIS (user ruling 2026-11-22 "yes", the k-fold of the forces and the stress; the forward's twin is
//! \c PairProjection::foldDomain).  Several blocks of a periodic density are one Gaussian lattice basis at several k: the functions'
//! values, gradients and pair products in real space are k-INDEPENDENT, and only each block's coefficient \f$D_{ij}e^{-ik\cdot R}\f$
//! differs.  A derivative that walks real space can therefore walk ONCE for all of them, if it knows which blocks share the walk and each
//! block's phase.  CLIENT: a per-block derivative face's folded variant (\c PointSetPullback::MakePointSetPullbackFolded), reached by an
//! abstract -> abstract cross-cast from the block it is handed.  Data-free.
class KFoldable
{
public:
    virtual ~KFoldable() = default;
    //! Blocks may fold together when their domains are EQUAL and NON-ZERO (the same lattice basis instance and grid ladder).
    virtual std::array<uint64_t,2> FoldDomain() const=0;
    //! My Bloch phase on the lattice translation \a n, \f$e^{ik\cdot R_n}\f$ in the convention my own collocation and gather use.
    virtual dcmplx FoldPhase(const ivec3_t& n) const=0;
};

//! One block of a fold: the block and its density matrix on the complex face (a real TRIM block's D widened; the caller keeps it alive).
struct KFoldTerm
{
    const KFoldable* block=nullptr;
    const chmat_t*   D    =nullptr;
};

} // namespace qchem::BasisSet

// File: BasisSet/DBCacheClient.C  The integral cache's KEY CONTRACT -- what a basis promises the cache.
//
// V1.20d: this struct is implemented by the PUBLIC face IrrepBasisSet_IDs, so it is public vocabulary and
// must not live in an .Internal. module (CLAUDE.md: never re-export an Internal module from a public one).
// The cache MECHANISM (IntegralsCache<T>, theCache) stays in qchem.BasisSet.Internal.DB_Cache -- only the
// contract is promoted.
module;
#include <cstddef>
#include <string>
export module qchem.BasisSet.DBCacheClient;

export namespace qchem::BasisSet {

//! A cache client supplies its own identity string.  The cache uses it verbatim as the per-basis key axis
//! and knows nothing about what is being cached (atoms, molecules, solids).  The contract: equal physics
//! (same exponents / angular momenta / contraction AND same centres / orientation) MUST give equal strings,
//! and any difference MUST give different strings -- the client owns ID assembly because it differs
//! completely for atoms, molecules and solids.  (See IrrepBasisSet_IDs for the default atom assembly and
//! PGData for the molecular, geometry-aware one.)
struct DBCacheClient
{
    virtual ~DBCacheClient() = default;
    virtual std::string BasisSetID() const = 0;
    //! The leading dimension this client expects its cached 2-centre matrices to have (= number of basis
    //! functions).  The cache cross-checks it on every hit/insert, so a BasisSetID() that is not specific
    //! enough (two differently-sized geometries colliding on one key) is caught right at the cache boundary
    //! -- with operator/ID/dims in hand -- instead of as a Cholesky segfault three layers down.
    //! Named CacheDim() (NOT size()) on purpose: size() collides with the VectorFunction/diamond hierarchy.
    //! Final overrider is the single bridge in IrrepBasisSet<T>.
    virtual size_t CacheDim() const = 0;
};

} //namespace

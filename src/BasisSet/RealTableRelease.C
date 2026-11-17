// File: BasisSet/RealTableRelease.C  A fit basis that can free its REAL blocks' tables for a phase that does not use them.
//
// WHY (OpenWork "Route B", lane H 2026-11-17).  Route A of the real-TRIM response runs the ground state on REAL TRIM blocks and the
// linear response on their complex SIBLINGS.  A δ fit caches one point-value table per orbital block, so during the response the
// ground state's real tables (TiO₂ k222: 1.44 GB) sat beside the siblings' complex ones (2.87 GB), and the response run's high-water
// rose above forceComplex's.  Nothing in the response asks for a real block's table: it is dead weight for the response's lifetime.
//
// THE PHASE IS EXPLICIT (doc/Pins.md pin 11): the client that knows a phase has started -- the facade, which builds the response --
// opens a \c RealTablesReleased scope and the fit frees the tables at once.  Never an automatic "nobody asked lately" eviction.
// A table is a pure function of fixed inputs, so one asked for again is REBUILT, the same.  Inside the phase such a rebuild undoes
// the saving, so it is counted, printed, and reported when the phase closes; after the phase it is the ordinary first-use fill.
module;
#include <cstddef>
export module qchem.BasisSet.RealTableRelease;

export namespace qchem::BasisSet
{

//! \brief The capability of a fit basis that caches per-REAL-block tables (a δ fit's \f$\Phi_{gi}=\chi_i(r_g)\f$ over a real TRIM block
//! of a complex run): free them for a phase that does not use them.  Reached by an abstract -> abstract cross-cast from the fit face;
//! use it through \c RealTablesReleased, never the two calls by hand.
class RealTableRelease
{
public:
    virtual ~RealTableRelease() = default;
    //! Free every real block's table (and whatever is built on it), and return the bytes freed.  Opens the phase.
    virtual size_t BeginRealTableRelease() const=0;
    //! Close the phase and return how many real tables were REBUILT inside it (0 = the release held for the whole phase).
    virtual size_t EndRealTableRelease() const=0;
};

//! \brief THE RELEASE PHASE, scoped: frees the real tables now, closes the phase at the end of the scope.  It writes what it did one
//! level into the caller's report context (the caller opens the scope that names the phase): \c realTableBytesReleased at the start,
//! \c realTablesRebuiltInPhase at the end.
class RealTablesReleased
{
public:
    explicit RealTablesReleased(const RealTableRelease& fit);
    ~RealTablesReleased();
    RealTablesReleased(const RealTablesReleased&)=delete;
    RealTablesReleased& operator=(const RealTablesReleased&)=delete;
    size_t BytesReleased() const;   //!< what the release freed
private:
    const RealTableRelease& itsFit;
    size_t                  itsBytes;
};

} // namespace

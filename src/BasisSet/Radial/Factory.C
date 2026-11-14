// File:: BasisSet/Factory.C  Interfaces for various basis set factories.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.BasisSet.Radial.Factory;
export import qchem.BasisSet;
export import qchem.ElectronConfiguration;


export namespace qchem::BasisSet::Radial
{
    enum class Type {Slater,Gaussian,BSpline6,BSpliner6,Gaussian_RKB,Slater_RKB};
    // Use these for custom basis sets.
    Real_BS* Factory(const nlohmann::json&,const ElectronConfiguration&);
    Real_BS* Factory(const nlohmann::json&,size_t Z); //Use atomic number Z to deduce the ElectronConfiguration
    // The object below is to define restricted collection of atom basis sets so that the global DB_Cache to not get too large
    // when running 100s of unit tests.  But these 'canned' versions are also generally usefull for production work.
    //
    // High:   As big and dense as you can get without becoming numercially unstable.
    // Medium: High with one smallest exponent and ~10-20% of the largest exponents removed. 
    //         Giving up core cusp but not exponent density
    // Low:    High with every other expoenent removed.  Giving up density but not core cusp.
    //
    // (The old N3/N5 tiny fixed even-tempered pools (3/5 exponents, ltrim=0) were REMOVED: they are ill-
    //  conditioned and invalid for SCF, and were repeatedly mistaken for solver/code bugs.  Low is the
    //  smallest SCF-grade pool -- use it for the "small deterministic basis" unit tests that used N3/N5.)
    enum class BasisSetAccuracy {Low,Medium,High};
    std::string BasisSetAccuracyStrs[]={"Low    ", "Medium ","High   "};
    Real_BS* Factory(BasisSetAccuracy, Type,size_t Z);
    //! As above, but the irreps come from an explicit \a ec instead of Atom_EC(Z): \a Z still sets the
    //! exponent pool (the element), while \a ec sets WHICH angular irreps the basis carries.  Needed when
    //! the occupied config differs from the neutral element -- e.g. a pseudo-ion (F- closed p^6 vs neutral
    //! F open p^5); for a neutral atom pass Atom_EC(Z) and the result is identical to the Z overload.
    Real_BS* Factory(BasisSetAccuracy, Type, size_t Z, const ElectronConfiguration& ec);

    //! \brief A radial basis that builds its Coulomb (Hartree) fit at a LADDER of richnesses, so a caller can size the fit to an
    //! accuracy target (the atom facade's policy, doc/Pins.md pin 31).  The rungs are NESTED: each one's functions include the rung
    //! below's, so the fit space only grows and the Dunlap energy error \f$-\tfrac12(\Delta\rho|\Delta\rho)\f$ only shrinks up the ladder.
    //!  - rung 0: \c CreateCDFitBasisSet's own fit, the first (s) block's exponents x 2;
    //!  - rung 1: every block's exponents x 2, merged into one s set (an l-shell density r^{2l} e^{-2 zeta r} peaks where the s block's
    //!    exponents are thin: U's 5f costs the rung-0 fit 32 mHa, lane B 2026-11-14);
    //!  - rung r >= 2: rung 1 with 2^{r-1} - 1 geometric points inserted between each pair of neighbours.
    //! Slater and Gaussian NR bases carry it; a basis without a Coulomb fit (B-spline, Dirac) does not.
    class CoulombFitLadder
    {
    public:
        virtual ~CoulombFitLadder() = default;
        //! The fit at \a rung, owned by the caller.  THROWS for \a rung outside [0, kCoulombFitRungs).
        virtual rFIT_CD_ABS* MakeCoulombFitRung(int rung) const=0;
        static constexpr int kCoulombFitRungs=6;
    };
}






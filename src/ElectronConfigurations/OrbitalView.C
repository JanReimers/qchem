// File: ElectronConfigurations/OrbitalView.C  The OccupationPolicy's view of one orbital block (DIP).
//
// V1.11 increment 3.  The policy (this library) needs to read orbital occupations, degeneracies and
// orthonormal-basis coefficients for its MOM reference capture and overlap scores -- but qcOrbitals sits
// ABOVE qcElConfig in the link DAG (qcOrbitals -> qcChargeDensity -> ... -> qcElConfig), so the policy
// cannot import it.  This is the CLAUDE.md dependency-inversion example verbatim: the LOWER library owns
// the abstract face, the higher library's class (TOrbitals) derives from it, and instances passed down
// arrive as the abstraction.  The policy never learns what an Orbitals is.
module;
#include <cstddef>
export module qchem.ElectronConfiguration.OrbitalView;
import qchem.Types;   // vec_t<T> (the orthonormal-basis coefficient column)

export namespace qchem {

//! \brief What the occupation policy may read from one orbital block: per-orbital occupation, capacity,
//! and the orthonormal-basis coefficient column (metric = I, so overlaps are plain conjugate dots).
//! Implemented by \c Orbitals::TOrbitals<T>; consumed by \c OccupationPolicy<T>.
template <class T> class OrbitalView
{
public:
    virtual ~OrbitalView() = default;
    virtual size_t   NumOrbitals()        const=0;
    virtual double   Occupation (size_t i) const=0;   //!< the whole level's occupation (g·f)
    virtual double   Degeneracy (size_t i) const=0;   //!< the level capacity g
    virtual vec_t<T> CoeffPrime (size_t i) const=0;   //!< orthonormal-basis coefficients C'
};

} // namespace qchem

// File: Structure/Lattice_3D/ReciprocalLattice.C  The reciprocal lattice dual to a Lattice_3D.
module;
#include <vector>
#include <iosfwd>
#include <cassert>
#include <typeinfo>   // std::bad_cast (thrown by the reference dynamic_cast in the pry-out helpers)
export module qchem.ReciprocalLattice;
export import qchem.UnitCell;
import qchem.Structure;   // the abstract base the pry-out helpers take
import qchem.Streamable;
import qchem.Math;        // FourPi (the Coulomb/Poisson kernel prefactor)
import qchem.Vector3D;    // rvec3_t + its dot product (G*G)

namespace qchem {

//! \brief The reciprocal lattice dual to a (direct) Lattice.
//!
//! Holds the reciprocal UnitCell — whose cell matrix is \f$B = 2\pi A^{-\top}\f$,
//! so \f$b_i\cdot a_j = 2\pi\delta_{ij}\f$ — and enumerates the reciprocal
//! lattice vectors \f$G = B m\f$ within a cutoff.  Symbol/units conventions are
//! documented in Lattice.C.  Obtain one from Lattice::Reciprocal().
export class ReciprocalLattice
    : public virtual Streamable
{
public:
    explicit ReciprocalLattice(const UnitCell& reciprocalCell) : itsCell(reciprocalCell) {}

    //! The reciprocal cell (its lattice matrix is \f$B\f$).
    const UnitCell& GetCell() const {return itsCell;}

    //! Reciprocal lattice vectors as integer index triples \f$m\f$ (\f$G = B m\f$)
    //! with \f$|G| \le\f$ Gmax.
    std::vector<ivec3_t> GetGVectors(double Gmax) const {return itsCell.CellsInSphere(Gmax);}

    //! Magnitude \f$|G|\f$ for an index triple \f$m\f$.  [1/length]
    double GetGLength(const ivec3_t& m) const {return itsCell.GetDistance(m);}

    //! \brief The diagonal Coulomb (Poisson) kernel \f$4\pi/|G|^2\f$ for index triple \a m (\f$G=B\,m\f$),
    //! with \f$m=0\to 0\f$ (the dropped neutralising background).  This is the reciprocal-space physics BOTH a
    //! plane-wave basis (its Hartree matrices) and a plane-wave density (the SAD seed's \f$V_H\f$) apply, so it
    //! lives here on the lattice metric \f$B\f$ -- not on the basis, which only borrows it.
    double CoulombKernel(const ivec3_t& m) const
    {
        if (m.x==0 && m.y==0 && m.z==0) return 0.0;
        rvec3_t G=itsCell.ToCartesian(rvec3_t(m));   // G = B m
        return FourPi/(G*G);
    }

    std::ostream& Write(std::ostream& os) const {return itsCell.Write(os);}

private:
    UnitCell itsCell; //!< Reciprocal cell: its lattice matrix is \f$B = 2\pi A^{-\top}\f$.
};

//---------------------------------------------------------------------------------
// Free pry-out helpers that downcast an abstract Structure to the concrete periodic UnitCell.
//
// Same pattern -- and the same rationale -- as Symmetry::Atom::Getl: the reciprocal lattice is wanted
// in several STRUCTURE-NEUTRAL libraries (the SCF iterator, density mixing, seeding), each of which was
// repeating its own dynamic_cast + assert.  These keep the cast AND its error handling in ONE place, and
// let the call site say what it WANTS ("give me the reciprocal lattice") instead of how to get it.
//
// The Get* forms THROW std::bad_cast on a non-periodic Structure (the reference dynamic_cast does it).
// isPeriodicCell() is the non-throwing probe, for a caller that must degrade gracefully rather than
// fail -- e.g. the density-mixer factory, which falls back to linear D-mixing.

//! Is this Structure a periodic cell (i.e. can the Get helpers below answer)?  Never throws.
export inline bool isPeriodicCell(const Structure* st)
{
    return dynamic_cast<const UnitCell*>(st)!=nullptr;
}

//! The periodic cell itself.  \throws std::bad_cast if \a st is not periodic.
export inline const UnitCell& GetUnitCell(const Structure* st)
{
    assert(st);
    return dynamic_cast<const UnitCell&>(*st);
}

//! The reciprocal cell \f$B=2\pi A^{-\top}\f$.  \throws std::bad_cast if \a st is not periodic.
export inline UnitCell GetReciprocalCell(const Structure* st)
{
    return GetUnitCell(st).MakeReciprocalCell();
}

//! The reciprocal lattice (the Poisson/Kerker metric).  \throws std::bad_cast if \a st is not periodic.
export inline ReciprocalLattice GetReciprocalLattice(const Structure* st)
{
    return ReciprocalLattice(GetReciprocalCell(st));
}

} // namespace qchem
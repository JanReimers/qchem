// File: Structure/Lattice_3D/ReciprocalLattice.C  The reciprocal lattice dual to a Lattice_3D.
module;
#include <vector>
#include <iosfwd>
export module qchem.ReciprocalLattice;
export import qchem.UnitCell;
import qchem.Streamable;

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

    std::ostream& Write(std::ostream& os) const {return itsCell.Write(os);}

private:
    UnitCell itsCell; //!< Reciprocal cell: its lattice matrix is \f$B = 2\pi A^{-\top}\f$.
};

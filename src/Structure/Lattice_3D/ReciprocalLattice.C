// File: Structure/Lattice_3D/ReciprocalLattice.C  The reciprocal lattice dual to a Lattice_3D.
module;
#include <vector>
#include <iosfwd>
export module qchem.ReciprocalLattice;
export import qchem.UnitCell;
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

} // namespace qchem
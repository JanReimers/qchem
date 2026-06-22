// File: Structure/UnitCell.C  Unit cell for a lattice.
module;
#include <iosfwd>
#include <vector>
export module Structure.UnitCell;
export import qchem.Types;
import qchem.Structure;
import qchem.Matrix3D;
import qchem.Streamable;

//! \brief One periodic unit cell: the cell geometry (lattice matrix \f$A\f$)
//! plus the atom basis it holds.  Symbol and units conventions (\f$A,a,\alpha,
//! M,B,G,k,\dots\f$) are documented in Lattice.C.
export class UnitCell
    : public virtual Structure
    , private Molecule //Holds the atom basis.
{
public:
    UnitCell(double a);                                   //!< Cubic, edge \f$a\f$ (a.u.).
    UnitCell(double a, double b, double c, double α, double β, double γ); //!< Lengths a.u., angles degrees.
    UnitCell(const Matrix3D<double>& A);                  //!< From a cell matrix \f$A\f$ (columns are the lattice vectors \f$a_i\f$).

    using Molecule::Insert;
    using Molecule::GetNumAtoms;

    UnitCell MakeReciprocalCell() const;                  //!< Reciprocal cell, \f$B = 2\pi A^{-\top}\f$.

    //! Cartesian (a.u.) position of a point given in fractional cell coordinates: \f$ r = A f \f$.
    rvec3_t ToCartesian(const rvec3_t& f) const;

    double      GetCellVolume     (                  ) const; //!< \f$|\det A|\f$ (a.u.\f$^3\f$).
    double      GetMinimumCellEdge(                  ) const; //!< \f$\min_i |a_i|\f$ (a.u.).
    double      GetDistance       (const rvec3_t& f  ) const; //!< \f$\sqrt{f^\top M f}\f$ for fractional \f$f\f$ (a.u.).
    vec3_t<int> GetNumCells       (double MaxDistance) const; //!< Cells per axis to cover a sphere of radius MaxDistance.

    //! Integer cell-index triples \f$n\f$ with \f$\lVert A n\rVert \le\f$ MaxDistance.  Pure
    //! cell geometry, so it serves both direct (R) and reciprocal (G) lattices.
    std::vector<vec3_t<int>> CellsInSphere(double MaxDistance) const;

    std::ostream&  Write(std::ostream&) const;

private:
    Matrix3D<double> itsA; //!< Cell matrix; columns are the lattice vectors \f$a_i\f$ (a.u.).
    Matrix3D<double> itsM; //!< Metric tensor \f$M = A^\top A\f$ (a.u.\f$^2\f$), cached.
};



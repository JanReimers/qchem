// File: Structure/UnitCell.C  Unit cell for a lattice.
module;
#include <iosfwd>
#include <vector>
#include <functional>
export module qchem.UnitCell;
export import qchem.Types;
import qchem.Structure;
import qchem.Matrix3D;
import qchem.Streamable;
import qchem.Mesh;        // qcMesh::Mesh / MeshParams (the CreateIntegrationMesh override)

namespace qchem {

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

    bool isFinite() const override {return false;}   //!< A periodic cell is NOT finite (Vnn -> Ewald).

    //! Form-factor sum PER CELL VOLUME: \f$\frac1\Omega\sum_a f(Z_a)\f$ -- the periodic G=0 background density
    //! (the finite-structure sum, normalised by \f$\Omega\f$).  Lets the PP G=0 alignment read the geometry
    //! without a Structure->UnitCell downcast or an LSP-violating CellVolume() on finite structures.
    double SumFormFactors(const std::function<double(int Z)>& f) const override
    {return Structure::SumFormFactors(f)/GetCellVolume();}

    //! A periodic cell's real-space integration mesh: a UNIFORM grid of \c mp.nUniform points per axis at
    //! cell-fractional midpoints (weight \f$\Omega/n^3\f$ each) -- the working lattice mesh for real-space PP
    //! quadrature.  (Plane-wave DFT integrates in G-space on the basis's own grid, so it never asks for this;
    //! an adaptive unit-cell Becke grid is a future refinement.)
    qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    //! \brief Add an atom of nuclear charge \a Z at FRACTIONAL cell coordinates \a f (\f$r=Af\f$).
    //! Convenience over Insert(new Atom(Z, ToCartesian(f))) so a crystal basis can be specified in
    //! cell coordinates (the natural way to give a diamond/FCC two-atom basis).
    void AddAtom(int Z, const rvec3_t& f);

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

    std::ostream&  Write(std::ostream&) const override;

private:
    Matrix3D<double> itsA; //!< Cell matrix; columns are the lattice vectors \f$a_i\f$ (a.u.).
    Matrix3D<double> itsM; //!< Metric tensor \f$M = A^\top A\f$ (a.u.\f$^2\f$), cached.
};

//! \brief Face-centred-cubic primitive cell of conventional cubic lattice constant \a a: the primitive
//! lattice vectors (columns of \f$A\f$) are the half-face diagonals \f$\tfrac a2(0,1,1),\,\tfrac a2(1,0,1),
//! \,\tfrac a2(1,1,0)\f$; cell volume \f$a^3/4\f$.  Add the atom basis with AddAtom in fractional
//! coordinates (e.g. diamond: \f$(0,0,0)\f$ and \f$(\tfrac14,\tfrac14,\tfrac14)\f$).
export class FCCUnitCell : public UnitCell
{
public:
    explicit FCCUnitCell(double a);
};



} // namespace qchem
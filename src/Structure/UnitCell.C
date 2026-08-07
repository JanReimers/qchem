// File: Structure/UnitCell.C  Unit cell for a lattice.
module;
#include <iosfwd>
#include <string>
#include <vector>
#include <functional>
#include <memory>
export module qchem.UnitCell;
export import qchem.Types;
import qchem.Structure;
import qchem.Matrix3D;
import qchem.Streamable;
import qchem.Mesh;        // qcMesh::Mesh / MeshParams (the CreateIntegrationMesh override)
import qchem.Symmetry.Lattice_3D.Fold;   // SymOp {W|τ} (the site-adapted Becke mesh's imposed ops)

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

    //! Polymorphic deep copy preserving periodicity (cell matrix + atom basis) -- so the facade keeps a
    //! UnitCell as a UnitCell.  (A FCCUnitCell clones to a plain UnitCell; it carries no state or behaviour
    //! beyond the cell geometry its base already holds, so the copy is geometrically identical.)
    std::shared_ptr<Structure> Clone() const override {return std::make_shared<UnitCell>(*this);}

    bool isFinite() const override {return false;}   //!< A periodic cell is NOT finite (IonIon -> Ewald).

    //! Form-factor sum PER CELL VOLUME: \f$\frac1\Omega\sum_a f(Z_a)\f$ -- the periodic G=0 background density
    //! (the finite-structure sum, normalised by \f$\Omega\f$).  Lets the PP G=0 alignment read the geometry
    //! without a Structure->UnitCell downcast or an LSP-violating CellVolume() on finite structures.
    double SumFormFactors(const std::function<double(int Z)>& f) const override
    {return Structure::SumFormFactors(f)/GetCellVolume();}

    //! A periodic cell's real-space integration mesh; \c mp.cellKind selects the quadrature.
    //! \c UnitCellKind::Uniform -- a UNIFORM grid of \c mp.nUniform points per axis at cell-fractional
    //! midpoints (weight \f$\Omega/n^3\f$ each): the working lattice mesh for real-space PP quadrature.
    //! (Plane-wave DFT integrates in G-space on the basis's own grid, so it never asks for this.)
    //! \c UnitCellKind::Becke -- the atom-centred periodic Becke fuzzy-Voronoi quadrature (dense radial
    //! near each nucleus, negligible-cost diffuse tails): the near-ideal grid for pointwise-nonlinear,
    //! sharp-at-the-core fields (XC).  The image series is \f$\varepsilon\f$-converged by magnitude
    //! screening; no radius parameter exists.  See doc/GPWPlan1.md "Becke XC grid".
    qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    //! \brief The same integration mesh, made INVARIANT under \a ops (doc/SymmetryUpgradePlan.md §6a) --
    //! the §3-imposed-symmetry sibling of the overload above.  Invariance is the T2 precondition for
    //! pointwise star-averaging and folding.
    //!
    //! HOW it is achieved is this class's business, not the caller's, and it follows \c mp.cellKind:
    //!  - \c Becke -- SITE-ADAPTED (W2b), invariant BY CONSTRUCTION: per atom ORBIT the representative
    //!    carries a site-group-invariant angular set (\c MakeInvariantAngularMesh under its Cartesian
    //!    stabilizer) and each symmetry partner the op-rotated copy.  No group-average growth pass.
    //!  - \c Uniform -- the generic mesh group-averaged onto the op orbits (W1's \c MakeInvariant; a
    //!    no-op dedup when the grid is already commensurate).
    //!
    //! It is NOT called \c CreateSiteAdaptedBeckeMesh any more (user, 2026-08-07): "Becke" merely repeated
    //! what \c mp.cellKind already says -- and the old body ASSERTED that value, i.e. the name was carrying
    //! a preconditon the parameter should carry -- while "SiteAdapted" named one of the two strategies,
    //! which is an implementation detail the caller should not have to choose.  \a ops present is what
    //! distinguishes this from the plain overload; that is what the signature should say.
    //! \note Still \c UnitCell-only, because \a ops is a \c Lattice_3D type.  A site-adapted MOLECULAR
    //! mesh is equally meaningful (user), and hoisting this to \c Structure needs a structure-neutral op
    //! type first -- see the worklist item.
    qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&,
                                       const std::vector<Symmetry::Lattice_3D::SymOp>& ops) const;

    //! \brief Add an atom of nuclear charge \a Z at FRACTIONAL cell coordinates \a f (\f$r=Af\f$).
    //! Convenience over Insert(new Atom(Z, ToCartesian(f))) so a crystal basis can be specified in
    //! cell coordinates (the natural way to give a diamond/FCC two-atom basis).  \a spinFlip stamps the
    //! atom's collinear magnetic-configuration flip bit (Atom::itsSpinFlip -- the spin-SAD seed's -m
    //! sublattice, doc/SCFSeedingPlan.md §10), so an AFM cell is specified atom by atom.
    void AddAtom(int Z, const rvec3_t& f, bool spinFlip=false);

    UnitCell MakeReciprocalCell() const;                  //!< Reciprocal cell, \f$B = 2\pi A^{-\top}\f$.

    //! The cell matrix \f$A\f$ (columns are the lattice vectors \f$a_i\f$, a.u.).  Exposed for the space-group
    //! detector (SpaceGroup::Detect wants \f$A\f$ + the fractional atom basis).
    const Matrix3D<double>& GetCellMatrix() const {return itsA;}

    //! Cartesian (a.u.) position of a point given in fractional cell coordinates: \f$ r = A f \f$.
    rvec3_t ToCartesian(const rvec3_t& f) const;

    //! Fractional cell coordinates of a Cartesian point: \f$ f = A^{-1} r \f$ (the inverse of \c ToCartesian).
    //! Used to place a real-space point on the FFT grid (its grid index is \f$N\,f\f$, modulo \f$N\f$).
    rvec3_t ToFractional(const rvec3_t& r) const;

    double      GetCellVolume     (                  ) const; //!< \f$|\det A|\f$ (a.u.\f$^3\f$).
    double      GetMinimumCellEdge(                  ) const; //!< \f$\min_i |a_i|\f$ (a.u.).
    double      GetMaximumCellEdge(                  ) const; //!< \f$\max_i |a_i|\f$ (a.u.); the axis that binds an isotropic uniform-mesh division count.
    double      GetDistance       (const rvec3_t& f  ) const; //!< \f$\sqrt{f^\top M f}\f$ for fractional \f$f\f$ (a.u.).
    vec3_t<int> GetNumCells       (double MaxDistance) const; //!< Cells per axis to cover a sphere of radius MaxDistance.

    //! Integer cell-index triples \f$n\f$ with \f$\lVert A n\rVert \le\f$ MaxDistance.  Pure
    //! cell geometry, so it serves both direct (R) and reciprocal (G) lattices.
    std::vector<vec3_t<int>> CellsInSphere(double MaxDistance) const;

    //! \brief Announce a UNIFORM \f$N_1\times N_2\times N_3\f$ grid over THIS cell: one console line
    //! (\c [uniform grid] <key>: N, points, \f$\Delta r_i=|a_i|/N_i\f$ a.u.) plus a \c grids.<key>
    //! run-report entry with the same fields (+ \c kind=Uniform, and \c eCut when \a eCut >= 0).
    //! THE single owner of uniform-cell-grid reporting: CreateIntegrationMesh's uniform branch AND the
    //! FFT raster engines (whose grids never pass through the mesh factory -- they recover this cell via
    //! MakeReciprocalCell) both delegate here, so every uniform grid announces in one format.  (The
    //! Reporting dependency stays in the implementation -- json in this interface leaks nlohmann
    //! reachability into every importer and collides with <ranges> in test TUs.)
    void EmitUniformGridReport(const std::string& key, const vec3_t<int>& N, double eCut=-1.0) const;

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
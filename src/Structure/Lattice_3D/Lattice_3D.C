// File: Structure/Lattice_3D/Lattice_3D.C Define a 3D infinite lattice.
module;
#include <vector>
#include <iosfwd>
#include <memory>

export module qchem.Lattice_3D;
export import qchem.Structure;
export import qchem.UnitCell;
export import qchem.ReciprocalLattice;
export import qchem.KMesh;
export import qchem.Symmetry.Lattice_3D.SpaceGroup;   // the crystal OWNS its space group (2026-08-01, plan §6a)

namespace qchem {

//! \brief A crystal lattice: a UnitCell repeated on a Bravais lattice, together
//! with the Brillouin-zone sampling used for plane-wave / band-structure work.
//!
//! A lattice is a block of unit cells.  A molecule-type structure inserted at
//! construction is repeated in every cell; the spatial extent then fixes the
//! number of k-points in the first zone and the plane-wave normalization.
//!
//! \par Symbol and units conventions (Hartree atomic units throughout)
//! Lengths are in Bohr \f$a_0\f$, energies in Hartree, reciprocal-space
//! quantities in \f$a_0^{-1}\f$.  Capitals denote \f$3\times3\f$ matrices (or
//! the lattice vectors that form their columns); lower case denotes lengths,
//! continuous points, or dimensionless fractional coordinates.
//!
//! \par Direct (real) space
//!  - \f$A\f$ — cell matrix; columns are the primitive lattice vectors
//!    \f$a_1,a_2,a_3\f$.  A Cartesian position is \f$r = A f\f$.  [length]
//!  - \f$a_i\f$ — primitive lattice vectors (columns of \f$A\f$).  [length]
//!  - \f$a,b,c\f$ — cell edge lengths \f$|a_i|\f$.  [length]
//!  - \f$\alpha,\beta,\gamma\f$ — cell angles, \f$\alpha=\angle(a_2,a_3)\f$,
//!    \f$\beta=\angle(a_1,a_3)\f$, \f$\gamma=\angle(a_1,a_2)\f$.  [rad, dimensionless]
//!  - \f$M\f$ — metric tensor \f$M = A^\top A\f$, so \f$\lVert r\rVert^2 = f^\top M f\f$.  [length\f$^2\f$]
//!  - \f$R\f$ — a direct lattice translation, \f$R = A n,\; n\in\mathbb{Z}^3\f$.  [length]
//!  - \f$r\f$ — a Cartesian position.  [length]
//!  - \f$f,n\f$ — fractional cell coordinates (\f$n\f$ integer); \f$r = A f\f$.  [dimensionless]
//!
//! \par Reciprocal space (the \f$2\pi\f$ lives in \f$B\f$, so \f$b_i\cdot a_j = 2\pi\delta_{ij}\f$)
//!  - \f$B\f$ — reciprocal cell matrix, \f$B = 2\pi A^{-\top}\f$; columns \f$b_i\f$.  [1/length]
//!  - \f$b_i\f$ — primitive reciprocal lattice vectors.  [1/length]
//!  - \f$G\f$ — a reciprocal lattice vector \f$G = B m,\; m\in\mathbb{Z}^3\f$; the
//!    plane-wave label in \f$e^{iG\cdot r}\f$.  [1/length]
//!  - \f$k\f$ — a crystal momentum, i.e. a continuous point in the Brillouin zone.  [1/length]
//!  - \f$K\f$ — a sampled \f$k\f$-point (a node of the Monkhorst–Pack mesh).  [1/length]
//!
//! \par Cutoffs
//! \f$E_{\max}\f$ (Hartree) bounds the plane-wave set \f$\{\,G : \tfrac12\lVert k+G\rVert^2 < E_{\max}\,\}\f$.
//
export class Lattice_3D
{
public:
    typedef std::shared_ptr<Structure> st_t;
    Lattice_3D(const UnitCell&, const Vector3D<int>&);                //Empty unit cell.
    
    const  UnitCell& GetUnitCell() const
    {
        return itsUnitCell;
    }

    //! \brief The crystal's atomic structure -- a copy of the unit cell's atom basis (Cartesian a.u.) as
    //! an abstract Structure.  The Hamiltonian's external term consumes this: the lattice carries all the
    //! information needed to build the external potential's structure factor.
    std::shared_ptr<const Structure> GetStructure() const;

    double GetLatticeVolume() const
    {
        return GetNumUnitCells()*itsUnitCell.GetCellVolume();
    }
    //! The dual reciprocal lattice (\f$B = 2\pi A^{-\top}\f$).  Apply the energy /
    //! \f$|G|\f$ cutoff later via ReciprocalLattice::GetGVectors().
    ReciprocalLattice Reciprocal() const;
    //! Monkhorst–Pack Brillouin-zone sampling using itsLimits as the divisions.  \a shift is the fractional
    //! grid-step offset (\f$0\f$ = Γ-centred; \f$½\f$ = the classic MP offset, i.e. CP2K's default for even grids).
    KMesh       MakeKMesh(const rvec3_t& shift={0,0,0}) const {return KMesh(itsLimits,shift);}
    ivec3_t     GetLimits() const {return itsLimits;}

    //! \brief The crystal's SPACE GROUP -- detected from the cell matrix + fractional atom basis on
    //! first call and cached (geometry is immutable, so the cache is shared across copies).  The
    //! UnitCell \f$\to\f$ {A, fractional sites} adapter lives HERE, common to EVERY lattice basis
    //! flavour (PW/GPW/APW/LAPW) instead of private to one factory (doc/SymmetryUpgradePlan.md §6a:
    //! qcStructure is where Mesh and Symmetry meet).  What a run IMPOSES from this group remains the
    //! §3 run-level policy, owned by the driver -- the lattice only answers what the symmetry IS.
    const Symmetry::Lattice_3D::SpaceGroup& GetSpaceGroup(double tol=1e-4) const;

    //! \brief The COLLINEAR magnetic (Shubnikov) group of this lattice under the per-atom decoration
    //! \a spins (Shubnikov S3, doc/SymmetryUpgradePlan.md §7 step 7): the same UnitCell \f$\to\f$
    //! {A, sites} adapter as \c GetSpaceGroup, with \c AtomSite::spin filled from \a spins (parallel to
    //! the cell's atom order -- assemble it with \c ChargeDensity::MagneticDecoration so it follows the
    //! seed's own species resolution).  Forwards to \c SpaceGroup::ShubnikovOps, which re-enumerates
    //! every \f$\tau\f$-coset (the anti-translation of a magnetically doubled cell -- the \f$m_1=-m_2\f$
    //! mirror -- is invisible to the detected one-coset-per-\f$W\f$ group).
    std::vector<Symmetry::Lattice_3D::SymOp> ShubnikovOps(const std::vector<int>& spins,
                                                          double tol=1e-4) const;

    size_t    GetNumSites     () const;
    size_t    GetNumBasisSites() const;
    size_t    GetNumUnitCells () const;

    size_t    GetSiteNumber   (const rvec3_t&  ) const;
    size_t    GetBasisNumber  (const rvec3_t&  ) const;
    size_t    GetBasisNumber  (size_t SiteNumber) const;

    Vector3D<int> GetCellCoord   (const rvec3_t&  ) const;
    rvec3_t         GetCoordinate  (size_t SiteNumber) const;
    void          SplitCoordinate(const rvec3_t& r, rvec3_t& basis, Vector3D<int>& cell) const;

    rvec_t      GetDistances    (size_t NumShells) const;
    rvec3vec_t  GetBonds        (size_t BasisNumber, double distance) const;
    rvec3vec_t  GetBondsInSphere(size_t BasisNumber, double distance) const;
    std::vector<ivec3_t>  GetCellsInSphere(double distance) const;

    std::ostream& Write(std::ostream&) const;

private:
    size_t      Find(const rvec3_t&               ) const; //Search unit cell.
    size_t      Find(double,const std::vector<double>&) const;
    rvec3vec_t  GetSuperCells(double MaxDistance) const;
    rvec3_t        GetBasisVector(size_t BasisNumber  ) const;

    UnitCell       itsUnitCell;  //Unit cell dimensions, no atoms.
    Vector3D<int>  itsLimits;    //Number of unit cell in each direction.
    double         itsTolerence; //Positions closer than this are considered the same.
    //! Lazily-detected space group (shared_ptr: copies of an immutable lattice share the cache).
    mutable std::shared_ptr<const Symmetry::Lattice_3D::SpaceGroup> itsSpaceGroup;
};


} // namespace qchem
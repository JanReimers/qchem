// File: BasisSet/Lattice_3D/PlaneWave_IBS.C  Plane-wave irrep basis set for one k-point.
//
// A complex (dcmplx) Orbital_1E_IBS whose functions are the normalised plane waves
// e^{i(k+G).r}/sqrt(V) for the reciprocal lattice vectors G in the cutoff set
// { G : 1/2 |k+G|^2 < Ecut }.  The wave-vector k labels the Bloch (translational)
// symmetry of the block (BlochFactory).
//
// Milestone 1 (empty lattice, V_ext=0): only Overlap (identity) and Kinetic
// (diagonal |k+G|^2) are exercised; the eigenvalues of 1/2*Kinetic reproduce the
// exact free-electron ladder 1/2 |k+G|^2.  See doc/PlaneWavePlan.md (sections 2.1, 3, 6).
module;
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

export module qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Internal.IrrepBasisSetImp;   // IrrepBasisSetImp<T>: GetSymmetry/GetSymt/GetIrrep
export import qchem.ReciprocalLattice;             // ctor takes a ReciprocalLattice (carries the B cell)
import qchem.Structure;
import qchem.Types;

export namespace BasisSet::Lattice_3D
{

//! \brief Plane-wave basis for a single k-point: the normalised waves
//! \f$ e^{i(k+G)\cdot r}/\sqrt V \f$ over the cutoff set \f$\{G:\tfrac12|k+G|^2<E_{cut}\}\f$.
class PlaneWave_IBS
    : public virtual BasisSet::Orbital_1E_IBS<dcmplx> // overlap/kinetic/nuclear + symmetry interface
    , public         BasisSet::IrrepBasisSetImp<dcmplx> // supplies GetSymmetry/GetSymt/GetIrrep + itsSymmetry
{
public:
    //! \param recip   the reciprocal lattice (its UnitCell matrix is \f$B=2\pi A^{-\top}\f$).
    //! \param N       Brillouin-zone grid divisions (context for the integer k-label).
    //! \param kIndex  integer k-label; the fractional crystal momentum is \f$k = kIndex/N\f$.
    //! \param Ecut    plane-wave energy cutoff (Hartree): keep \f$G\f$ with \f$\tfrac12|k+G|^2<E_{cut}\f$.
    PlaneWave_IBS(const ReciprocalLattice& recip, const ivec3_t& N,
                  const ivec3_t& kIndex, double Ecut);

    virtual size_t GetNumFunctions() const {return itsG.size();}

    //! Reciprocal-index label \f$m\f$ of basis function \a i (its plane wave is \f$e^{i(k+G)\cdot r}\f$,
    //! \f$G = B\,m\f$).  This integer triple is the defining quantum number of the plane wave.
    ivec3_t GetGIndex(size_t i) const {return itsG[i];}

    // 1E integral building blocks (no 1/2 on Kinetic -- the Hamiltonian applies it).
    virtual csmat_t MakeOverlap () const;                  //!< Identity (PWs orthonormal over the cell).
    virtual csmat_t MakeKinetic () const;                  //!< Diagonal \f$|k+G|^2 = \langle p^2\rangle\f$.
    virtual csmat_t MakeNuclear (const Structure*) const;  //!< Zero stub (see TODO -- milestone 2.3).

    //! \brief Assemble \f$ \langle G|V|G'\rangle = \tilde V(G-G') \f$ from a caller-supplied G-space
    //! potential, keyed by the reciprocal-index difference \f$\Delta m = m(G)-m(G')\f$.
    //!
    //! This is the reusable G-space potential assembly: the milestone-2.2 separable cosine and the
    //! milestone-2.3 nuclear structure factor are both just particular \f$\tilde V\f$ suppliers.
    //! \note Returns smat_t (symmetric) -- correct for a real, even \f$\tilde V\f$ (the cosine).  A
    //! genuinely Hermitian-but-not-symmetric \f$\tilde V\f$ (complex structure factor) will need hmat_t.
    csmat_t MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    // VectorFunction: the plane-wave values / gradients at a point.
    virtual cvec_t     operator() (const rvec3_t& r) const;
    virtual cvec3vec_t Gradient   (const rvec3_t& r) const;

    virtual std::string Name      () const {return "PlaneWave";}
    virtual std::string BasisSetID() const; // geometry-aware cache key (N, kIndex, Ecut, nG)

    virtual std::ostream& Write(std::ostream&) const;

private:
    rvec3_t GetGCartesian(const ivec3_t& m) const; //!< \f$ G = B\,m \f$ in Cartesian a.u.

    ReciprocalLattice    itsRecip;  //!< Reciprocal cell (matrix \f$B\f$); source of G and |k+G|.
    ivec3_t              itsN;      //!< BZ grid divisions (for the cache key).
    ivec3_t              itsKIndex; //!< Integer k-label (for the cache key).
    rvec3_t              itsK;      //!< Fractional crystal momentum \f$k=kIndex/N\f$.
    double               itsEcut;   //!< Energy cutoff (Hartree).
    double               itsVolume; //!< Direct cell volume \f$V\f$ (for the \f$1/\sqrt V\f$ norm).
    std::vector<ivec3_t> itsG;      //!< Surviving reciprocal index triples \f$m\f$ (\f$G=Bm\f$).
};

} //namespace

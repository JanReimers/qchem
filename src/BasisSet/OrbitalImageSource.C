// File: BasisSet/OrbitalImageSource.C  "Under this crystal op, my orbitals become these": the capability face the
// little-group machinery needs from a Bloch block (doc/KAxisPlan.md P4 + §6 ruling: the basis ANSWERS, it does not
// hand out its shells, centres or Cartesian powers for a caller to rotate).
//
// Consumers: D5-a (the response's k on the little group of q: rotate a block's δD onto its star), D5-b / CK-3a (rebuild
// an unstored block C_{gk} = M C_k from an IBZ state), and the P4 gates.  A block that is a Bloch sum over atom-centred
// shells answers through Symmetry::Lattice_3D::BuildBlochOperationRep; a plane-wave block would answer with a G
// permutation and phase.
module;
export module qchem.BasisSet.OrbitalImageSource;
export import qchem.Symmetry.SymOp;   // SymOp {W|tau,sigma}, fractional direct frame
export import qchem.Types;            // mat_t<dcmplx>

export namespace qchem::BasisSet
{
class OrbitalImageSource
{
public:
    virtual ~OrbitalImageSource() = default;
    //! \brief The matrix \f$M\f$ with \f$\hat g\,\phi_{k,\mu}=\sum_\nu M(\nu,\mu)\,\phi_{k',\nu}\f$: this block's orbitals
    //! carried by \a g onto the block at \f$k'=W^{-T}k\f$ (\c Symmetry::Lattice_3D::ImageOfK; any stored block equal to it
    //! mod G -- the Bloch functions are periodic in k).  Orbital coefficients map as \f$C_{k'}=MC_k\f$; an operator matrix
    //! as \f$F_k=M^\dagger F_{k'}M\f$ when \a g is a symmetry of it.  \a g's spin action is NOT applied: it swaps the spin
    //! BLOCKS, which the caller owns.  THROWS when \a g is not a symmetry of this basis on its cell.  Complex on every
    //! block (a real TRIM block's image is real only when \f$k'\f$ is TRIM too, which the caller can check).
    virtual mat_t<dcmplx> OrbitalImage(const Symmetry::SymOp& g) const = 0;
};
}

// File: BasisSet/Lattice_3D/LAPW_IBS.C  Linearized Augmented Plane Wave (LAPW) basis for one k-point.
//
// Lineage B, second IBS (see doc/PlaneWavePlan.md).  LAPW fixes APW's energy-dependent (nonlinear)
// secular equation: inside the muffin-tin sphere the single energy-dependent radial function is replaced
// by a fixed linear combination of u_l(r) = u_l(r,E_l) and its energy derivative udot_l = du_l/dE at a
// FIXED linearization energy E_l.  Matching BOTH value and radial derivative at R fixes the two
// coefficients (a_l, b_l), and H and O become energy-INDEPENDENT -> a single generalized eigenproblem
// H c = e O c yields a whole band at once (a true band solver, unlike APW's per-energy determinant).
//
// Empty lattice (V=0): u_l(r)=j_l(q r), q=sqrt(2 E_l), and udot_l(r)=(r/q) j_l'(q r).  H is built from the
// symmetric-gradient form 1/2 int grad.phi* grad.phi' (manifestly Hermitian; avoids the u/udot Hamiltonian
// asymmetry), O from int phi* phi'.  The eigenvalues reproduce the free-electron ladder 1/2|k+G|^2 up to
// the LAPW linearization error (smallest near E_l).
//
// Like APW_IBS, LAPW_IBS shares the IrrepBasisSet<dcmplx> base with PlaneWave_IBS but is its own class
// with its own (generalized-eigenproblem) interface -- lineage B is a sequence of such IBS types.
module;
#include <iosfwd>
#include <string>
#include <vector>

export module qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.BasisSet.IrrepBasisSet;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
export import qchem.ReciprocalLattice;
import qchem.Types;

export namespace BasisSet::Lattice_3D
{

//! \brief Linearized Augmented Plane Wave basis for one k-point (single origin muffin-tin sphere).
class LAPW_IBS
    : public virtual BasisSet::IrrepBasisSet<dcmplx>
    , public         BasisSet::IrrepBasisSetImp<dcmplx>
{
public:
    //! \param Ecut interstitial plane-wave cutoff; \param Rmt muffin-tin radius; \param lmax angular
    //! cutoff; \param Elin linearization energy E_l (>0) for the radial functions u_l, udot_l.
    LAPW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
             double Ecut, double Rmt, size_t lmax, double Elin);

    virtual size_t GetNumFunctions() const {return itsG.size();}

    //! \brief Energy-independent overlap O(k) and Hamiltonian H(k) (empty lattice, V=0).  Both Hermitian
    //! (real for the single origin sphere).  Solve the generalized problem H c = e O c for the band.
    chmat_t MakeOverlap()     const;
    chmat_t MakeHamiltonian() const;

    // VectorFunction: interstitial plane-wave value (in-sphere augmentation not evaluated pointwise).
    virtual cvec_t     operator() (const rvec3_t& r) const;
    virtual cvec3vec_t Gradient   (const rvec3_t& r) const;

    virtual std::string Name      () const {return "LAPW";}
    virtual std::string BasisSetID() const;
    virtual std::ostream& Write(std::ostream&) const;

private:
    // Assemble both H and O in one radial-quadrature + matching pass (they share a_l,b_l).
    void Assemble(chmat_t& H, chmat_t& O) const;

    ReciprocalLattice    itsRecip;
    rvec3_t              itsK;
    double               itsVolume;
    double               itsRmt;
    size_t               itsLmax;
    double               itsElin;
    std::vector<ivec3_t> itsG;
};

} //namespace

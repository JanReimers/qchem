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
// Because H and O are energy-independent at fixed E_l, LAPW_IBS honours the SAME Orbital_1E_IBS<dcmplx>
// interface as atoms / molecules / plane waves: MakeOverlap, MakeKinetic (<p^2>, no 1/2), MakeNuclear.
// The Hamiltonian is then assembled the usual way, H = 1/2 Kinetic + Nuclear, and solved as the
// generalized problem H c = e O c.  (The energy-dependent friction was specific to *self-consistent*
// APW; fixing E_l removes it -- the augmented basis functions become fixed objects.)
module;
#include <iosfwd>
#include <string>
#include <vector>

export module qchem.BasisSet.Lattice_3D.LAPW_IBS;
import qchem.BasisSet.Orbital_1E_IBS;             // Orbital_1E_IBS<dcmplx>: Overlap/Kinetic/Nuclear
import qchem.BasisSet.Internal.IrrepBasisSetImp;
export import qchem.ReciprocalLattice;
import qchem.Structure;
import qchem.Types;

export namespace BasisSet::Lattice_3D
{

//! \brief Linearized Augmented Plane Wave basis for one k-point (single origin muffin-tin sphere).
class LAPW_IBS
    : public virtual BasisSet::Orbital_1E_IBS<dcmplx>
    , public         BasisSet::IrrepBasisSetImp<dcmplx>
{
public:
    //! \param Ecut interstitial plane-wave cutoff; \param Rmt muffin-tin radius; \param lmax angular
    //! cutoff; \param Elin linearization energy E_l for the radial functions u_l, udot_l; \param Z
    //! nuclear charge of the (single, origin) atom -- muffin-tin potential V=-Z/r inside the sphere,
    //! V=0 in the interstitial.  Z=0 is the empty lattice (free radial functions, analytic Bessel).
    //! In an augmented method the potential defines the basis, so Z is fixed at construction.
    LAPW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
             double Ecut, double Rmt, size_t lmax, double Elin, double Z=0.0);

    virtual size_t GetNumFunctions() const {return itsG.size();}

    // Orbital_1E_IBS building blocks (energy-independent at fixed E_l; Hermitian, real for one origin
    // sphere).  Assemble H = 1/2 Kinetic + Nuclear and solve the generalized problem H c = e O c.
    virtual chmat_t MakeOverlap() const;                 //!< <phi_i|phi_j>
    virtual chmat_t MakeKinetic() const;                 //!< <p^2> = <phi_i|-grad^2|phi_j> (NO 1/2)
    virtual chmat_t MakeNuclear(const Structure*) const; //!< <phi_i|V|phi_j>, V=-Z/r (uses the ctor Z)

    // VectorFunction: interstitial plane-wave value (in-sphere augmentation not evaluated pointwise).
    virtual cvec_t     operator() (const rvec3_t& r) const;
    virtual cvec3vec_t Gradient   (const rvec3_t& r) const;

    virtual std::string Name      () const {return "LAPW";}
    virtual std::string BasisSetID() const;
    virtual std::ostream& Write(std::ostream&) const;

private:
    // Assemble the kinetic <p^2> block, the overlap, and the nuclear <V> block in one radial-quadrature
    // + matching pass (they share the radial tables and the matching coefficients a_l,b_l).
    void Assemble(chmat_t& Kp2, chmat_t& O, chmat_t& Vnuc) const;

    ReciprocalLattice    itsRecip;
    rvec3_t              itsk;
    double               itsVolume;
    double               itsRmt;
    size_t               itsLmax;
    double               itsElin;
    double               itsZ;       //!< Nuclear charge (muffin-tin V=-Z/r); 0 = empty lattice.
    std::vector<ivec3_t> itsG;
};

} //namespace

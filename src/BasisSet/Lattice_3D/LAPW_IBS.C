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
    //! cutoff; \param Elin linearization energy E_l for the radial functions u_l, udot_l; \param Znuc
    //! nuclear charge of the (single, origin) atom -- muffin-tin potential V=-Znuc/r inside the sphere,
    //! V=0 in the interstitial.  Znuc=0 is the empty lattice (free radial functions, analytic Bessel).
    //! In an augmented method the potential defines the basis, so Znuc is fixed at construction.
    LAPW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
             double Ecut, double Rmt, size_t lmax, double Elin, double Znuc=0.0);

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
    //! One angular channel l, expressed in the two-function radial basis {u_l, udot_l}.  The three
    //! operator blocks are 2x2 symmetric matrices in that basis, so a muffin-tin matrix element is the
    //! quadratic form  c_i . (block . c_j)  in the matching coefficients (see CombineBlocks).
    struct RadialBlock
    {
        rmat2d_t boundary;   //!< [[u, udot],[u', udot']] at r=Rmt; det = Wronskian, and it solves the matching
        rmat2d_t overlap;    //!< <u_a|u_b>
        rmat2d_t kinetic;    //!< <p^2> = <-grad^2>   (NO 1/2)
        rmat2d_t potential;  //!< <u_a|V|u_b>,  V=-Z/r
    };
    //! The augmentation of one plane wave k+G: its value+slope matching coefficients (a_l,b_l) per l.
    //! (The k+G geometry itself lives in Internal::KPlusG.)
    struct AugmentedWave
    {
        std::vector<rvec2d_t> c;      //!< matching coefficients (a_l, b_l), one 2-vector per l
        explicit AugmentedWave(std::vector<rvec2d_t> c_) : c(std::move(c_)) {}
    };

    // The constructor runs these three acts once and stores the assembled blocks (below).
    std::vector<RadialBlock>   MuffinTinRadialBlocks() const;                              // Act 1
    std::vector<AugmentedWave> MatchAugmentation(const std::vector<RadialBlock>&) const;   // Act 2
    void                       CombineBlocks(const std::vector<RadialBlock>&,
                                             const std::vector<AugmentedWave>&);           // Act 3

    ReciprocalLattice    itsRecip;
    rvec3_t              itsk;
    double               itsVolume;
    double               itsRmt;
    size_t               itsLmax;
    double               itsElin;
    double               itsZnuc;    //!< Nuclear charge (muffin-tin V=-Znuc/r); 0 = empty lattice.
    std::vector<ivec3_t> itsG;
    chmat_t              itsOvlp, itsKp2, itsVnuc;   //!< assembled in the constructor (Act 3)
};

} //namespace

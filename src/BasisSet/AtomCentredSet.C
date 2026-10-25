// File: BasisSet/AtomCentredSet.C  An ATOM-CENTRED FUNCTION SET, and the overlap service an orbital basis
// offers against one -- one module, because neither half means anything without the other (the
// Orbital_PP_IBS pattern).  doc/HubbardProjectorPlan.md D1 + D2.
//
// WHY.  A +U projector used to be built from the site's OWN orbital shells (a contraction over them), so the
// projector was whatever those shells could hold, and it changed with the basis (β-MnO₂: the O 2p projector
// captured 0.976 of the pseudo-atom orbital; one diffuse shell moved the gap 0.128 -> 0.237 eV).  QE defines
// it as an atomic function, independent of the basis.  The precedent is in this tree: a separable PP's
// projectors β are NOT basis functions; a Gaussian basis forms b_i = <χ_i|β> analytically.  A +U projector
// is the same kind of object, minus the weights: the 2l+1 functions R(r) Y_lm of one manifold on one site.
//
// It is per SITE (a manifold is), not per species like \c SpeciesProjectorSet; no common base was forced:
// the two share the radial convention (Math::Gaussian with \f$r^l\f$ carried by the owner) and nothing else.
//
// THE ANGULAR CONVENTION: real solid harmonics, \c Math::SphericalShell(l)[l+m] for m = -l..l in that order
// (l=1: y, z, x; l=2: xy, yz, 3z²-r², xz, x²-y²), each scaled by a POSITIVE factor to unit norm on the sphere
// (the GPW KB projectors' RealYlm table up to a sign per m).  A consumer that needs the site-group action on the
// set (+U slots) builds it from (centre, l) and this convention: a SphericalShellRep over the same raw harmonic.
module;
#include <cstddef>
#include <stdexcept>
#include <vector>
export module qchem.BasisSet.AtomCentredSet;
export import qchem.Math.Gaussian;            // Math::Gaussian -- one c r^{2n} e^{-alpha r^2} term
export import qchem.BasisSet.BareCoulombSource; // the set's own bare integrals (ACBN0 over a pseudo-atomic manifold, HubbardProjectorPlan D6)
import qchem.Types;                           // rvec3_t, mat_t<T>

export namespace qchem::BasisSet
{
//! \brief A set of atom-centred functions, grouped into MANIFOLDS: manifold \f$M\f$ is the \f$2l_M+1\f$ functions
//! \f$\chi_{M,m}(r)=R_M(|r-C_M|)\,Y_{l_Mm}(\widehat{r-C_M})\f$, \f$m=-l..l\f$, each unit-normalised
//! (\f$\int R_M^2r^2dr=1\f$).  The functions of one manifold are mutually orthonormal; functions of different
//! manifolds are whatever their geometry makes them.  Columns of any matrix over the set run manifold by
//! manifold, \f$m\f$ ascending.
class AtomCentredSet
{
public:
    virtual ~AtomCentredSet() = default;
    virtual size_t  Count () const=0;            //!< the number of manifolds
    virtual rvec3_t Centre(size_t M) const=0;    //!< \f$C_M\f$ (Cartesian, bohr)
    virtual int     L     (size_t M) const=0;    //!< \f$l_M\f$
    size_t Size(size_t M) const {return size_t(2*L(M)+1);}
    //! \f$\sum_M(2l_M+1)\f$: the column count of a matrix over the set.
    size_t NumFunctions() const {size_t n=0; for (size_t M=0;M<Count();M++) n+=Size(M); return n;}
    //! The first column of manifold \a M.
    size_t Offset(size_t M) const {size_t n=0; for (size_t K=0;K<M;K++) n+=Size(K); return n;}
};

//! \brief CAPABILITY: each radial is a closed sum of Gaussians, \f$R_M(r)=\sum_tc_t\,r^{\,l_M+2n_t}e^{-\alpha_tr^2}\f$
//! (the \f$r^l\f$ carried by the owner, as \c SpeciesProjectorSet_Gaussian) -- so a Gaussian basis forms its
//! overlaps with the set analytically.  Reached by the abstract->abstract cross-cast.
class AtomCentredSet_Gaussian : public virtual AtomCentredSet
{
public:
    virtual ~AtomCentredSet_Gaussian() = default;
    virtual const std::vector<Math::Gaussian>& AsGaussians(size_t M) const=0;
};

//! \brief The value type: manifolds added one by one, each radial NORMALISED on entry (so a caller hands over
//! a shape, e.g. an atomic orbital's raw contraction, and the set's contract holds by construction).
//! It answers its OWN bare two-electron integrals (\c BareCoulombSource, over set column indices) for
//! functions on ONE centre -- what ACBN0 needs over a manifold (HubbardProjectorPlan D6, increment 4).
class GaussianAtomCentredSet
    : public virtual AtomCentredSet_Gaussian
    , public virtual BareCoulombSource
{
public:
    //! Add manifold (\a centre, \a l) with radial shape \a radial (any overall scale; \f$r^l\f$ implied).
    //! Returns its index.  Throws on an empty or zero-norm radial.
    size_t Add(const rvec3_t& centre, int l, std::vector<Math::Gaussian> radial);
    //! \f$\int_0^\infty R(r)^2r^2dr\f$ for \f$R=\sum_tc_tr^{l+2n_t}e^{-\alpha_tr^2}\f$ (exact).
    static double RadialNorm2(int l, const std::vector<Math::Gaussian>& radial);

    virtual size_t  Count ()         const override {return itsM.size();}
    virtual rvec3_t Centre(size_t M) const override {return itsM[M].centre;}
    virtual int     L     (size_t M) const override {return itsM[M].l;}
    virtual const std::vector<Math::Gaussian>& AsGaussians(size_t M) const override {return itsM[M].radial;}
    //! \brief \f$(ab|cd)\f$ over set columns \a cols, all on ONE centre (throws otherwise: a two-centre integral
    //! has no business here).  EXACT, no quadrature: the one-centre multipole expansion
    //! \f$(ab|cd)=\sum_k\frac{4\pi}{2k+1}R^k(ab;cd)\sum_qG_{ab}^{kq}G_{cd}^{kq}\f$, with the real Gaunt factors
    //! \f$G\f$ from sphere-monomial integrals of the set's own harmonics and the Slater integrals \f$R^k\f$ of
    //! Gaussian radial densities in closed form (an incomplete beta \f$B_x(a,b)\f$ with INTEGER \f$b\f$, a finite
    //! sum).  Central cell, no lattice sum -- the \c BareCoulombSource contract.
    virtual ERI4Block BareCoulomb(const std::vector<size_t>& cols) const override;
private:
    struct Manifold { rvec3_t centre; int l; std::vector<Math::Gaussian> radial; };
    std::vector<Manifold> itsM;
};

//! \brief THE SERVICE (D2): "the overlap of every basis function with each function of an atom-centred set".
//! \f$B_{i,(M,m)}=\langle\chi_i|\chi_{M,m}\rangle\f$, \f$n\times\f$ \c NumFunctions(), in the block scalar \tparam T.
//! A periodic block answers for its Bloch functions against the set's Bloch sums at the block's \f$k\f$,
//! \f$B_{i,(M,m)}=\sum_Re^{ik\cdot R}\langle\chi_i|\chi_{M,m}(\cdot-R)\rangle\f$ -- Hermitian-consistent with its
//! \c MakeOverlap: if \f$\chi_{M,m}=\sum_jc_j\chi_j\f$ then \f$B=Sc\f$.  An integral a basis can compute knowing
//! nothing about Hubbard.  A block that cannot answer does not derive the face (the consumer throws).
//! Uncached (\c Make): the consumer holds what it builds from it.
template <class T> class AtomCentredOverlapSource
{
public:
    virtual ~AtomCentredOverlapSource() = default;
    virtual mat_t<T> MakeAtomCentredOverlap(const AtomCentredSet&) const=0;
    //! \brief The set's OWN overlap \f$O\f$ at this block's \f$k\f$ (D3): \f$O_{ab}=\sum_Re^{ik\cdot R}\langle\chi_a|\chi_b(\cdot-R)\rangle\f$,
    //! the functions' TRUE mutual overlap, lattice-summed (what QE's ortho-atomic \f$O^{-1/2}\f$ is of).  None of
    //! this block's functions enter: the block is asked only because the k convention and the lattice sums are its
    //! side's.  Compare \f$B^\dagger S^{-1}B\f$, the part of \f$O\f$ the block REPRESENTS.
    virtual hmat_t<T> MakeSetOverlap(const AtomCentredSet&) const=0;
};

} //namespace

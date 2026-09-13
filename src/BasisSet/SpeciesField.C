// File: BasisSet/SpeciesField.C  The ARGUMENT vocabulary of the species-field integral service.
//
// V1.2 (doc/CleanupCandidates.md): a pseudopotential asks a basis for two NEW TYPES OF INTEGRAL -- the matrix
// of a species-attached local radial field, and the matrix of a weighted set of species-attached radial x Y_lm
// projectors.  Both are integrals a basis can compute knowing nothing about pseudopotentials, so the ARGUMENT
// TYPES are defined here, structure-neutrally, and the service that takes them is Orbital_PP_IBS<T>.  The
// pseudopotential library then IMPLEMENTS these faces (its HGH/GTH models are species radial fields and
// projector sets) and DEPENDS ON this library -- the old edge, qcBasisSet -> qcPseudopotential, is gone.
//
// Every face here is data-free and abstract; the optional ones (_Gaussian, _R) are CAPABILITIES reached by
// the sanctioned abstract->abstract cross-cast, so a field that is not a closed sum of Gaussians simply does
// not carry that face and the consumer keeps its quadrature route.  Each capability DERIVES (virtually) from
// its core face -- "a projector set that also answers in r" -- so one reference to the capability serves a
// consumer that needs both, and the diamonds this makes are the house style.
module;
#include <vector>
export module qchem.BasisSet.SpeciesField;
import qchem.Math.Gaussian;   // Math::Gaussian -- one c r^{2n} e^{-alpha r^2} term

export namespace qchem::BasisSet
{

//! \brief A RANGE SPLIT of a species field: the CALLER's routing choice, not a property of the field.  A
//! periodic KS assembly folds a field's long (softened-Coulomb) tail into the Hartree Poisson and keeps only
//! the compact short remainder in the external term (CP2K's local-PP split, doc/GPWPlan.md 0e-PP); a
//! consumer with no such routing asks for Full.  A field with no natural split answers Full for the whole
//! and 0 for Short.
enum class FieldRange { Full, Long, Short };

//! \brief A radial scalar field attached to each SPECIES \f$Z\f$: \f$v_Z(r)\f$.  The argument of the local
//! integral \f$\langle i|\sum_a v_{Z_a}(|r-R_a|)|j\rangle\f$ -- the basis supplies the lattice sum, the
//! structure factor and the cell volume; the field supplies only its one-species radial shape.
//!
//! DUAL-SPECTRAL BY DESIGN: a field answers in \f$r\f$ (a real-space quadrature route) AND in \f$q\f$ (the
//! analytic-FT route, \f$\tilde v(q^2)=\int d^3r\,v(r)e^{-iq\cdot r}\f$, which a periodic basis needs for
//! box independence -- not a speed hack).  The closed-form models have both views; a tabulated one would
//! transform once.
class SpeciesRadialField
{
public:
    virtual ~SpeciesRadialField() {}
    //! \f$v_Z(r)\f$ in a.u. for the given range.
    virtual double ValueR(int Z, double r,  FieldRange) const=0;
    //! \f$\tilde v_Z(q^2)\f$, \f$q>0\f$, for the given range. [energy x volume]
    virtual double ValueQ(int Z, double q2, FieldRange) const=0;
    //! \brief The finite \f$q\to0\f$ limit of the NON-Coulomb content of the given range, \f$\int d^3r\,
    //! [v_Z(r)+Z_{ion}/r]\f$ -- the uniform shift a periodic lattice sum drops with its \f$G=0\f$ term and has
    //! to align to (the "\f$\alpha\f$" of a plane-wave total energy).  0 when the range has none.
    virtual double CellMeanQ(int Z, FieldRange) const=0;
};

//! \brief CAPABILITY: the field is a finite closed sum of Gaussians, \f$v_Z(r)=\sum_t c_t r^{2n_t}
//! e^{-\alpha_t r^2}\f$ for the given range (exact, not fitted).  A consumer holding a Gaussian basis can then
//! assemble \f$\langle\chi_i|v|\chi_j\rangle\f$ analytically (a 3-centre Gaussian overlap) and size its grids
//! from the sharpest \f$\alpha\f$.  Optional: reached by cross-cast from \c SpeciesRadialField.
class SpeciesRadialField_Gaussian : public virtual SpeciesRadialField
{
public:
    virtual ~SpeciesRadialField_Gaussian() {}
    virtual std::vector<Math::Gaussian> AsGaussians(int Z, FieldRange) const=0;
};

//! \brief A WEIGHTED SET of species-attached radial \f$\times\,Y_{lm}\f$ projectors, \f$\sum_p w_p
//! |\beta_p Y_{l_p m}\rangle\langle\beta_p Y_{l_p m}|\f$ per species: the argument of the separable
//! integral \f$\langle i|\sum_a\sum_p w_p|\beta^a_p\rangle\langle\beta^a_p|\,|j\rangle\f$ (Kleinman-Bylander
//! form).  \c Weight is already a per-projector SCALAR -- a model with a coupling matrix diagonalises it
//! before it gets here, so no matrix of couplings ever reaches a basis.  The \f$m\f$ sum is the basis's.
class SpeciesProjectorSet
{
public:
    virtual ~SpeciesProjectorSet() {}
    virtual size_t Count  (int Z)           const=0;   //!< projectors for species Z
    virtual int    L      (int Z, size_t p) const=0;   //!< angular momentum of projector p
    virtual double Weight (int Z, size_t p) const=0;   //!< its weight \f$w_p\f$ [energy]
    //! The reciprocal radial shape \f$\tilde\beta_p(q)\f$, \f$q=|k+G|\f$ (the view a plane-wave basis consumes).
    virtual double RadialQ(int Z, size_t p, double q) const=0;
};

//! \brief CAPABILITY: the real-space radial shape \f$\beta_p(r)\f$ (the view a real-space quadrature
//! consumes).  Its spherical-Bessel transform reproduces \c RadialQ: \f$\int_0^\infty\beta_p(r)j_l(qr)r^2dr
//! =\tilde\beta_p(q)/\sqrt{4\pi}\f$.
class SpeciesProjectorSet_R : public virtual SpeciesProjectorSet
{
public:
    virtual ~SpeciesProjectorSet_R() {}
    virtual double RadialR(int Z, size_t p, double r) const=0;
};

//! \brief CAPABILITY: \f$\beta_p(r)\f$ as a closed sum of Gaussians, \f$\beta_p(r)=\sum_t c_t r^{\,l+2n_t}
//! e^{-\alpha_t r^2}\f$ with \f$l=\f$ \c L(Z,p) carried by the owner -- so a Gaussian basis forms
//! \f$\langle\chi_i|\beta_pY_{lm}\rangle\f$ analytically.
class SpeciesProjectorSet_Gaussian : public virtual SpeciesProjectorSet
{
public:
    virtual ~SpeciesProjectorSet_Gaussian() {}
    virtual std::vector<Math::Gaussian> AsGaussians(int Z, size_t p) const=0;
};

} //namespace

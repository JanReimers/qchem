// File: BasisSet/Orbital_PP_IBS.C  The species-field integral service of an orbital basis, AND its argument
// vocabulary -- one module, because neither half means anything without the other (absorbed SpeciesField.C,
// user 2026-09-15).
//
// V1.2 (doc/CleanupCandidates.md): a pseudopotential asks a basis for two NEW TYPES OF INTEGRAL -- the matrix
// of a species-attached local radial field, and the matrix of a weighted set of species-attached radial x Y_lm
// projectors.  Both are integrals a basis can compute knowing nothing about pseudopotentials, so the ARGUMENT
// TYPES are defined here, structure-neutrally, beside the service that takes them (Orbital_PP_IBS<T>).  The
// pseudopotential library then IMPLEMENTS these faces (its HGH/GTH models are species radial fields and
// projector sets) and DEPENDS ON this library -- the old edge, qcBasisSet -> qcPseudopotential, is gone.
//
// Every argument face here is data-free and abstract; the optional ones (_Gaussian, _R) are CAPABILITIES reached
// by the sanctioned abstract->abstract cross-cast, so a field that is not a closed sum of Gaussians simply does
// not carry that face and the consumer keeps its quadrature route.  Each capability DERIVES (virtually) from
// its core face -- "a projector set that also answers in r" -- so one reference to the capability serves a
// consumer that needs both, and the diamonds this makes are the house style.
module;
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.BasisSet.Orbital_PP_IBS;
import qchem.Math.Gaussian;                   // Math::Gaussian -- one c r^{2n} e^{-alpha r^2} term (the _Gaussian faces answer in it)
import qchem.Structure;                       // Structure -- the positions and species the lattice sums run over
import qchem.Types;                           // hmat_t<T>
export import qchem.AD;                       // NoGradient: a family with no PP geometry derivative refuses (GeometryOptimizationPlan §13.2)
export import qchem.Outcome;

export namespace qchem::BasisSet
{
//------------------------------------------------------------------------------------------------------------
// THE ARGUMENT VOCABULARY.
//------------------------------------------------------------------------------------------------------------
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
    //! \name THE q-SLOPE beside \c ValueQ (doc/Records/CellHalf_Design.md §5; user ruling §6.2, 2026-11-10, the derivative version beside an
    //! existing method): \f$d\tilde v_Z/d(q^2)\f$ for the given range, \f$q>0\f$.  CLIENT: a periodic basis's strain derivative of the field's
    //! matrix (a homogeneous strain moves every \f$|\Delta G|^2\f$).  \c CanSlopeQ is false exactly when \c SlopeQ throws: a field with no
    //! closed-form slope keeps its matrix, and its stress is a named \c MissingVirials() entry.
    //!@{
    virtual bool   CanSlopeQ() const {return false;}
    virtual double SlopeQ(int, double, FieldRange) const
    {throw std::logic_error("SpeciesRadialField::SlopeQ: this field gives no q-slope (CanSlopeQ() says so before a stress is asked)");}
    //!@}
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

//! \brief A spherical Gaussian charge \f$q(\alpha/\pi)^{3/2}e^{-\alpha r^2}\f$ of attracting sign: its
//! potential on an electron is \f$-q\,\mathrm{erf}(\sqrt\alpha\,r)/r\f$.
struct GaussianCharge
{
    double q;       //!< the charge
    double alpha;   //!< its exponent
};

//! \brief CAPABILITY: the field's LONG range is the potential of a Gaussian charge, exactly (not fitted).  A
//! Gaussian basis then forms \f$\langle\chi_i|v_{long}|\chi_j\rangle\f$ as a 3-centre Coulomb integral.
//! Optional: reached by cross-cast from \c SpeciesRadialField.
class SpeciesRadialField_GaussianCharge : public virtual SpeciesRadialField
{
public:
    virtual ~SpeciesRadialField_GaussianCharge() {}
    virtual GaussianCharge AsGaussianCharge(int Z) const=0;
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
    //! \name THE SLOPE beside \c RadialQ (CellHalf_Design §5, as \c SpeciesRadialField::SlopeQ): \f$d\tilde\beta_p/dq\f$.  \c CanSlopeQ is
    //! false exactly when \c RadialQSlope throws.
    //!@{
    virtual bool   CanSlopeQ() const {return false;}
    virtual double RadialQSlope(int, size_t, double) const
    {throw std::logic_error("SpeciesProjectorSet::RadialQSlope: this projector set gives no q-slope (CanSlopeQ() says so before a stress is asked)");}
    //!@}
};

//! \brief CAPABILITY: the real-space radial shape \f$\beta_p(r)\f$ (the view a real-space quadrature
//! consumes).  Its spherical-Bessel transform reproduces \c RadialQ: \f$\int_0^\infty\beta_p(r)j_l(qr)r^2dr
//! =\tilde\beta_p(q)/\sqrt{4\pi}\f$.
class SpeciesProjectorSet_R : public virtual SpeciesProjectorSet
{
public:
    virtual ~SpeciesProjectorSet_R() {}
    virtual double RadialR(int Z, size_t p, double r) const=0;
    //! \brief The projector's SHARPNESS: the Gaussian exponent \f$\alpha\f$ whose \f$e^{-\alpha r^2}\f$ matches
    //! \f$\beta_p\f$'s sharpest radial feature (an HGH projector of radius \f$r_l\f$ answers \f$1/2r_l^2\f$).
    //! A real-space view is incomplete without its scale: a consumer that QUADRATURES \f$\langle\chi|\beta_p\rangle\f$
    //! on a uniform grid sizes that grid from \f$\alpha_{\max}+\alpha_\beta\f$, the integrand's exponent -- the
    //! same role \c qcMesh::XCMeshSharpness::alphaPP plays for the local field (doc/CleanupCandidates.md V2.5).
    virtual double SharpnessR(int Z, size_t p) const=0;
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

//------------------------------------------------------------------------------------------------------------
// THE SERVICE.
//------------------------------------------------------------------------------------------------------------
//! \brief The two integrals a species-attached external field asks of an orbital basis (V1.2; formerly
//! \c Pseudopotential::Integrals_Pseudo, whose arguments were pseudopotential types).  \tparam T the matrix
//! scalar.  A periodic basis assembles \f$\frac1\Omega\sum_a f(Z_a)\,e^{-i\Delta G\cdot\tau_a}\f$ from the
//! field's \f$q\f$ view; a real-space basis quadratures the \f$r\f$ view -- the SAME face, the basis owning HOW.
//! Nothing here is pseudopotential-specific: the arguments are a radial field and a projector set.
template <class T> class Orbital_PP_IBS
{
public:
    virtual ~Orbital_PP_IBS() {}
    //! Can this basis assemble \a field?  A term keeps its quadrature route when the answer is no: a closed-form
    //! Gaussian route needs a field that carries the closed-form faces, and a decorator asks what it decorates.
    virtual bool CanAssembleField     (const SpeciesRadialField&)  const {return true;}
    //! Can this basis assemble \a projectors?  As \c CanAssembleField.
    virtual bool CanAssembleProjectors(const SpeciesProjectorSet&) const {return true;}
    //! \f$\langle i|\sum_a v_{Z_a}(|r-R_a|)|j\rangle\f$ for the given RANGE of \a field (Hermitian).  The
    //! \f$\Delta G=0\f$ term of a periodic assembly is the dropped uniform shift (\c CellMeanQ aligns it).
    virtual hmat_t<T> MakeSpeciesFieldMatrix(const Structure*, const SpeciesRadialField& field, FieldRange) const=0;
    //! \f$\langle i|\sum_a\sum_p w_p|\beta^a_p\rangle\langle\beta^a_p|\,|j\rangle\f$ (Hermitian), the \f$m\f$
    //! sum done by the basis (a plane-wave basis via the \f$(2l+1)P_l(\cos\gamma)\f$ addition theorem).
    virtual hmat_t<T> MakeProjectorMatrix(const Structure*, const SpeciesProjectorSet&) const=0;
    //! The PULLBACKS of the two assemblies (GeometryOptimizationPlan §13.2): given the block's density matrix \a D (Hermitian),
    //! F[A] = Σ_ij Re(D_ij* ∂M_ij/∂R_A) for every atom A, the basis centres AND the species' sites moving.  CLIENT: the PP terms' forces.
    //! Each \c Can is false exactly when its \c Make refuses (NoGradient: a family or a field with no derivative yet).
    virtual bool CanPullbackSpeciesField(const SpeciesRadialField&)  const=0;
    virtual bool CanPullbackProjectors  (const SpeciesProjectorSet&) const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeSpeciesFieldPullback(const Structure*, const SpeciesRadialField& field, FieldRange,
                                                                        const hmat_t<T>& D) const=0;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeProjectorPullback(const Structure*, const SpeciesProjectorSet&, const hmat_t<T>& D) const=0;
    //! THE STRAIN DERIVATIVES beside them (doc/Records/CellHalf_Design.md §5; user ruling §6.2, 2026-11-10): \f$V_{ab}=\partial\,\mathrm{Re}
    //! \,\mathrm{Tr}(D^\dagger M)/\partial\varepsilon_{ab}\f$ under the homogeneous strain \f$R\to(1+\varepsilon)R\f$, \f$A\to(1+\varepsilon)A\f$ at
    //! fixed \a D, a GRADIENT in Ha, channel A (the centres) and channel B (the cell's own appearance: \f$\Omega\f$, \f$G\f$) together.
    //! CLIENT: the PP terms' virials.  Default: a named refusal (a molecule has no cell; a family with no strain derivative yet), and each
    //! \c Can is false exactly when its \c Make refuses.
    virtual bool CanVirialSpeciesField(const SpeciesRadialField&)  const {return false;}
    virtual bool CanVirialProjectors  (const SpeciesProjectorSet&) const {return false;}
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeSpeciesFieldVirial(const Structure*, const SpeciesRadialField&, FieldRange, const hmat_t<T>&) const
    {return Outcome<rmat3d_t,AD::NoGradient>::Fail({"this basis gives no species-field virial"});}
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeProjectorVirial(const Structure*, const SpeciesProjectorSet&, const hmat_t<T>&) const
    {return Outcome<rmat3d_t,AD::NoGradient>::Fail({"this basis gives no projector virial"});}
    //! \brief The PER-ANGULAR-CHANNEL decomposition of \c MakeProjectorMatrix: one matrix per \f$l\f$, summing
    //! exactly to the whole.  DIAGNOSTIC (doc/SphericalLatticePlan.md I0: the MnO ordering campaign reads
    //! \f$E^{(l)}=\mathrm{Tr}(DV^{(l)})\f$ per channel).  Default: the whole lumped under \f$l=-1\f$, so a basis
    //! that has not implemented the split still answers correctly (callers treat \f$-1\f$ as "unresolved").
    virtual std::map<int,hmat_t<T>> MakeProjectorMatrixByL(const Structure* st, const SpeciesProjectorSet& ps) const
    {
        std::map<int,hmat_t<T>> byL;
        byL.emplace(-1, MakeProjectorMatrix(st, ps));
        return byL;
    }
};

} //namespace

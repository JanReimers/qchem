// File: BasisSet/DeltaFit_IBS.C  The delta (identity) v_xc fit basis over a real-space mesh.
//
// The OTHER representation a periodic Vxc fit basis can take (PlaneWaveFit_IBS being the {G} one): the
// delta basis, whose n_pts functions are the mesh's own delta functions.  Everything it knows is the
// quadrature -- points and weights -- because that IS what a delta basis is: W_g[h] = w_g delta_gh, a
// family of weight vectors with nothing left over.  (It still HOLDS the whole FitQuadrature bundle, since
// that is what its factory builds and hands it, but nothing it answers reads the fold.)
//
// It exists so that CreateVxcFitBasisSet can make BOTH decisions -- which representation AND which points
// -- and return ONE object (doc/OpenWork.md, "separation of concerns in the XC terms").  Before this the
// mesh came back through a SECOND factory (CreateXCQuadrature) as a bare struct, and the Hamiltonian had
// to re-decide which term class to build in order to know which of the two answers to use; the fit basis
// and the grid were separate returns of one question.
//
// It ANSWERS, it does not hand out (user ruling 2026-08-22).  The mesh is private: what leaves this class
// is an integral over its own functions or a projection onto them -- arrays indexed by FUNCTION, in its own
// order.  No consumer knows whether these functions sit on an atom-centred Becke build or a uniform cell
// grid, which is exactly the fit/grid separation the design item is about.  Since 2026-08-22 NOTHING leaves:
// the mesh getter went with the Quadrature face, and since 2026-08-24 the two operations that were reaching
// PAST that -- the symmetrizations -- are gone too (the fold travels to their consumer with the mesh).
//
// AND SINCE 2026-08-23 NOTHING SAYS "POINT" EITHER (user).  A delta basis is a family of FUNCTIONS, so it
// presents exactly the interface a Gaussian auxiliary basis presents: GetNumFunctions() for the count,
// Overlap(f) = <delta_g|f> for the projection, Integrals() = <delta_g|1> for the per-function integrals.
// The old NumPoints/Sample/Integrate trio described the QUADRATURE behind them, and read as correct only
// because here n_functions == n_points, so the wrong accessor returned the right number.
module;
#include <array>
#include <functional>
#include <memory>
#include <mutex>
#include <cassert>
#include <map>
#include <ostream>
#include <string>
#include <type_traits> // std::is_same_v (the same-scalar table cache)
#include <utility>
#include <vector>
export module qchem.BasisSet.DeltaFit_IBS;
export import qchem.BasisSet.BasisCache;   // a block's constructor takes its parent's cache (IntegralCacheOwnershipPlan S3a)
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Fit_Types;   // FitQuadrature / VxcFit -- the fit-factory vocabulary                    // cFIT_SF_ABS + Integrals_Overlap3C + FitQuadrature
export import qchem.BasisSet.PointImages;          // PhiRealization / PhiChooser / BlochImageSource: the shared-Φ realization (2026-11-19)
export import qchem.BasisSet.Transition_DFT_IBS;   // Transition_Overlap3C: the (k+q, k) pair on my Φ tables (LinearResponsePlan B2)
export import qchem.BasisSet.GradOverlap3C;        // Integrals_GradOverlap3C: the GGA's ∇(χ_iχ_j) at my points (GGA stage 2)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.Symmetry;                                   // sym_t (the Bloch irrep)
import qchem.Symmetry.Group;   // Symmetry::GroupPtr / TrivialGroup (Forces, Group face step 3a)
import qchem.Types;                                      // dcmplx, rvec3_t, vec_t
import qchem.Mesh.Quadrature;                            // qcMesh::SiteIntegrals (the atomic-partition observable)
import qchem.SlotMap;                                    // InstanceSlotMap: the per-block caches, keyed by the block's instance (S3c)
import qchem.Reporting;                                  // EmitFold -- this basis announces its own star-average

export namespace qchem::BasisSet
{

//! \brief THE PERIODIC FIT AXIS'S TWO EXTRA FACES, which only a \f$\delta\f$ basis on the complex (Bloch) axis has
//! (2026-10-31): a REAL TRIM block contracted in real arithmetic (3c-3, \c Integrals_Overlap3C<double,dcmplx>) and the
//! \f$(k+q,k)\f$ transition pair (\c Transition_Overlap3C).  A base that EXISTS only for \c dcmplx, so the real
//! (molecular) instantiation does not carry capabilities it has no use for.  CRTP: the bodies are the δ basis's own
//! (its tables, its weights), reached through \a Self.
template <class Self, class T> class DeltaPeriodicFaces {};
template <class Self> class DeltaPeriodicFaces<Self,dcmplx>
    : public virtual Integrals_Overlap3C<double,dcmplx>  // a REAL TRIM block (3c-3), in real arithmetic
    , public virtual Transition_Overlap3C            // a (k+q, k) transition density on my points (the Becke/δ XC route)
{
public:
    //! \copydoc BasisSet::Integrals_Overlap3C::Overlap3C
    const PairProjection<double>& Overlap3C(const Orbital_DFT_IBS<double,dcmplx>& orb) const override {return Self_().RealTensor(orb);}
    //! \copydoc BasisSet::Transition_Overlap3C::TransitionForward
    //! Two Φ tables from the same per-block cache the ground state fills, so a response adds no table build.
    cvec_t TransitionForward(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                             const mat_t<dcmplx>& dD) const override {return Self_().TransitionForwardImp(bra,ket,dD);}
    //! \copydoc BasisSet::Transition_Overlap3C::TransitionAdjoint
    mat_t<dcmplx> TransitionAdjoint(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                    const cvec_t& v) const override {return Self_().TransitionAdjointImp(bra,ket,v);}
private:
    const Self& Self_() const {return static_cast<const Self&>(*this);}
};

//! \brief THE REAL (MOLECULAR) FIT AXIS'S EXTRA FACE (GGA stage 2, the molecular point route): the gradient 3-centre overlap
//! \f$\nabla(\chi_i\chi_j)(r_g)\f$ (\c Integrals_GradOverlap3C), which a GGA's \f$\nabla\rho\f$ and its gradient term of H go through.  A base
//! that EXISTS only for \c double in stage 2: the Bloch axis gains it in stage 3, on the shared-Φ images, so until then a complex δ basis
//! does not carry a capability it cannot serve (the gradient is a type: no face, no call).  CRTP: the bodies are the δ basis's own.
template <class Self, class T> class DeltaRealFaces {};
template <class Self> class DeltaRealFaces<Self,double>
    : public virtual Integrals_GradOverlap3C<double,double>
{
public:
    //! \copydoc BasisSet::Integrals_GradOverlap3C::GradOverlap3C
    const PairProjection<double>& GradOverlap3C(const Orbital_DFT_IBS<double,double>& orb) const override {return Self_().GradTensor(orb);}
private:
    const Self& Self_() const {return static_cast<const Self&>(*this);}
};

//! \brief The \f$\delta\f$-basis realization of the \f$v_{xc}\f$ fit face: it holds the finished XC quadrature and
//! answers nothing else.
//!
//! ONE CLASS TEMPLATE ON THE FIT AXIS'S SCALAR \a T (2026-10-31, user: "replace the Vxc fit basis with delta functions ...
//! if our OOD design (LSP) is working correctly, this should just work").  \a T is the scalar of the orbital blocks I am
//! contracted against -- \c dcmplx on the periodic (Bloch) path, \c double on the molecular one -- NOT what my own functions
//! are made of: a \f$\delta\f$ function and its weight are real either way, and widen to \a T in one place (\c AsT).
//! Until then the class was \c cFIT_SF_ABS only, so the real path THREW for VxcFit::Delta: the scalar axis had been welded
//! to the class (its old comment blamed the Bloch irrep, which is incidental -- every IrrepBasisSet names one).
//!
//! Symmetry: the CALLER's -- \f$\Gamma\f$ Bloch irrep on the periodic path (the mesh is cell-periodic and carries no crystal
//! momentum of its own, exactly as \c PlaneWaveFit_IBS's fit block is built at \f$k=0\f$), the molecular fit lineage's
//! trivial irrep on the molecular one (\f$v_{xc}\f$ is totally symmetric).
//! \c BasisSetID folds in the point count and the fold's orbit count -- the only things that distinguish two of these --
//! but nothing keys an integral cache on a \f$\delta\f$ basis today: it has no 3-centre tensor, because the quadrature
//! route never forms one.
template <class T> class DeltaFit_IBS
    : public virtual FIT_SF_ABS<T>                   // the neutral scalar-fit face (+ Integrals_Overlap3C<T,T>)
    , public         DeltaPeriodicFaces<DeltaFit_IBS<T>,T>   // dcmplx only: the real TRIM block + the transition pair
    , public         DeltaRealFaces<DeltaFit_IBS<T>,T>       // double only (stage 2): the GGA's gradient 3-centre overlap
    , public         IrrepBasisSetImp<T>             // GetSymmetry/GetSymt/GetIrrep
    , public         Symmetry::IdentityRepresentation<T>   // identity at E only: the points PERMUTE under an op, which is the
                                                            // orbit fold's business (a vector op), not a matrix image (S3)
{
    friend class DeltaPeriodicFaces<DeltaFit_IBS<T>,T>;
    friend class DeltaRealFaces<DeltaFit_IBS<T>,T>;
public:
    //! ORTHOGONAL, with \f$\langle\delta_g|\delta_g\rangle=w_g\f$ -- diagonal, so no metric SOLVE, which
    //! is what \c isOrtho asks.  (Not orthoNORMAL: a fit through this basis divides by that diagonal,
    //! \f$c_g=\langle\delta_g|f\rangle/w_g=f(r_g)\f$.)  Answered HERE since 2026-08-23 -- it was the last
    //! thing holding up a \f$\delta\f$-only FACE, and "my metric is diagonal" is a fact about this
    //! implementation, not a type other code needs to name.
    bool isOrtho() const override {return true;}
    //! \f$\delta\f$ functions at points impose nothing of their own: the POINTS carry any group (\c qcMesh::ImposedBy).
    qcMesh::ImposedSymmetry Imposes() const override {return qcMesh::ImposedSymmetry::None;}
    Symmetry::GroupPtr      Forces () const override {return Symmetry::TrivialGroup();}   //!< {e}: the points carry any group
    //! \brief How a REAL orbital block is evaluated at this basis's points, when its OWN basis set's factory says so (R1.0 step 2,
    //! user ruling 2026-11-14).  The atom's blocks are radial: the factory (\c Radial::BasisSet_HF) hands a closure that fills
    //! \f$\Phi_{gi}=R_i(|r_g|)\f$ from each block's radial face, the m count living in the occupation, so the atom SCF never asks a
    //! block's 3-D face -- which is display work, to become an honest \f$Y_{lm}\f$ per (i, m).  Unset (every other lineage): the block's
    //! own \c op(r).  Concrete wiring in the factory, never a capability cast here.
    using RealPointValues=std::function<mat_t<double>(const Orbital_1E_IBS<double>&, const rvec3vec_t&)>;
    //! \brief Where a block's k-independent lattice IMAGES come from, when its basis set's factory says so (the shared \f$\Phi\f$, lane S
    //! 2026-11-19): the periodic factory hands a closure answering a block's \c BlochImageSource face, or null for a block without one.
    //! Concrete wiring in the factory, never a capability cast here.  Empty (every non-periodic lineage): no images, dense tables.
    using ImageSourceOf=std::function<const BlochImageSource*(const IrrepBasisSet_IDs&)>;
    //! \a q is the FINISHED quadrature (mesh + fold + Shubnikov tags) from \c MakeXCQuadrature, and \a sym the irrep the caller's lineage gives a fit block (see the class note).
    //! \a realValues: see \c RealPointValues (empty = each block's \c op(r)).  \a imageSourceOf: see \c ImageSourceOf.
    DeltaFit_IBS(const std::shared_ptr<BasisCache>& cache, FitQuadrature q, const sym_t& sym, RealPointValues realValues={},
                 ImageSourceOf imageSourceOf={})
        : IrrepBasisSetImp<T>(sym, cache)
        , itsQuad(std::move(q))
        , itsRealValues(std::move(realValues))
        , itsImageSourceOf(std::move(imageSourceOf))
    {
        assert(itsQuad.GetMesh() && itsQuad.GetMesh()->size()>0 && "DeltaFit_IBS: a delta basis IS its mesh -- it cannot be empty");
        // The invariants of the bundle, checked HERE, once, where it becomes an object -- they used to be
        // re-checked by the consumer that received the loose struct.
        // The mesh<->fold consistency asserts that used to sit here are GONE: qcMesh::FoldedMesh's
        // constructor checks them once, for every consumer, and THROWS rather than asserting -- an assert
        // is compiled out under NDEBUG, which is where production lives (2026-09-08).
        // Providers self-report (user ruling 2026-08-16): the star-average this basis will apply to rho on
        // EVERY iteration announces itself at birth, armed or not.  nOps is known only for a magnetic fold
        // (sigmas run parallel to the ops); 0 = "not reported", NOT "not armed".
        qchem::report::EmitFold("XC mesh", itsQuad.NumSpinOps(), itsQuad.GetMesh()->size(),
                                itsQuad.NumOrbits(),
                                (itsQuad.NumSpinOps()==0) ? std::string() : std::string("magnetic (Shubnikov)"));
    }

    //! \copydoc BasisSet::Integrals_Overlap3C::Overlap3C
    //! \f$\langle\chi_i|\delta_g|\chi_j\rangle=w_g\overline{\chi_i(r_g)}\chi_j(r_g)\f$, built from the
    //! orbital basis's own \c op(r) -- no orbital-side 3-centre machinery is asked for or needed.  The SAME-scalar block
    //! here; on the periodic axis a real TRIM block takes \c DeltaPeriodicFaces' overload, through the same templated
    //! body (its \f$\Phi\f$ table is typed by the ORBITAL block, so a real block contracts in real arithmetic).
    const PairProjection<T>& Overlap3C(const Orbital_DFT_IBS<T,T>& orb) const override;

    // ---- the integrals over my own functions, one entry per FUNCTION (FIT_SF_ABS) ---------------------
    //! \copydoc BasisSet::FIT_SF_ABS::OverlapDiagonal
    //! \f$\langle\delta_g|\delta_g\rangle=w_g\f$: orthogonal, NOT orthonormal -- so a general fit
    //! through this basis divides by these, giving \f$c_g=w_g f_g/w_g=f(r_g)\f$, the point values.
    vec_t<T> OverlapDiagonal() const override {return AsT(itsQuad.GetMesh()->Weights());}
    //! \copydoc BasisSet::FIT_SF_ABS::Charge
    //! \f$\langle\delta_g|1\rangle=\int\delta_g\,d^3r=w_g\f$.  Same numbers as \c OverlapDiagonal here
    //! and NOT the same question: that one is \f$\langle\delta_g|\delta_g\rangle\f$ (a metric), this is
    //! \f$\langle\delta_g|1\rangle\f$ (an integral).  They coincide because \f$\delta\f$ is idempotent
    //! under this quadrature -- a fact about this representation, not a shared definition.
    vec_t<T> Charge() const override {return AsT(itsQuad.GetMesh()->Weights());}
    //! \copydoc BasisSet::FIT_SF_ABS::Overlap
    //! \f$\langle\delta_g|f\rangle=w_g f(r_g)\f$ -- I evaluate the field at MY points (the field's own
    //! bulk fast path does the work) and weight it.  The caller sees only a vector indexed by FUNCTION.
    //! The field samples itself ON my mesh (\c SampleOn): a Fourier-series field then asks the mesh HOW (a NUFFT on a Becke
    //! build, the exact inverse FFT on the uniform raster), and every other field is evaluated at the points, as before.
    vec_t<T> Overlap(const ScalarFunction<double>& f) const override
    {
        const rvec_t  v=f.SampleOn(*itsQuad.GetMesh());
        const rvec_t& w=itsQuad.GetMesh()->Weights();
        assert(v.size()==w.size());
        vec_t<T> p(v.size());
        for (size_t g=0; g<v.size(); g++) p[g]=T(w[g]*v[g]);
        return p;
    }

    // NO Symmetrize / SymmetrizeSpin since 2026-08-24 (user): the star-average was never a fit-basis
    // question -- this class was merely the only object holding the orbit fold.  The fold and the
    // Shubnikov tags now travel to the XC strategy WITH the mesh, through the same factory out-parameter,
    // and it applies the free SymmetrizeValues / SymmetrizeValuesSigned itself.  Same cure, and the same
    // reason, as SiteIntegrals leaving here in increment 5.

    //! \brief THE FACTORY'S SETUP STEP for the shared \f$\Phi\f$ (not a face: the factory holds the concrete object it just built).  Counts
    //! both realizations' memory for the orbital blocks \a cblocks (complex) and \a rblocks (real TRIM) -- the dense tables exactly, the image
    //! records from an exact count on a subsample -- asks \a choose (the policy layer's door) which realization to take, and on \c Images builds
    //! the shared records now, at setup, so the run's console reads in run order.  On \c Dense nothing is built here: the tables fill as they
    //! always did (lazily, per block, bit-identical).  Returns the realization taken.
    PhiRealization ChoosePhiRealization(const std::vector<const Orbital_1E_IBS<dcmplx>*>& cblocks,
                                        const std::vector<const Orbital_1E_IBS<double>*>&  rblocks, const PhiChooser& choose);

    //! One \f$\delta\f$ function per mesh point -- ~100k of them, so anything reasoning about a fit basis
    //! by COUNTING functions (grid sizing, memory reports, cache dimensions) must not choke on that.
    size_t GetNumFunctions() const override {return itsQuad.GetMesh()->size();}

    // NO op(r)/Gradient, and no stub either: since 2026-08-22 IrrepBasisSet does not promise pointwise
    // evaluation (doc/CleanupCandidates.md R1.0 step 1), so a basis whose functions are DISTRIBUTIONS
    // simply does not offer it.  This class used to carry two throwing overrides for exactly that reason
    // -- the interface demanded an answer it could not honestly give.  Nothing asks: the delta route's
    // fit coefficients ARE the field's values at the mesh points, and it has no expansion to plot, so it
    // does not implement FieldEvaluator either.

    std::string   Name      () const override {return "DeltaFit";}
    std::string   BasisSetID() const override
        {return Name()+"|npts="+std::to_string(itsQuad.GetMesh()->size())
                      +"|orbits="+std::to_string(itsQuad.NumOrbits());}
    std::ostream& Write     (std::ostream& os) const override
        {return os << Name() << " fit IBS: " << itsQuad.GetMesh()->size() << " delta functions on the XC mesh";}

private:
    //! The face is templated on the scalar of the ORBITAL blocks I contract against (Bloch => dcmplx), not
    //! on what my own functions are made of -- and a \f$\delta\f$ function and its weight are real.  So
    //! every integral over my own functions widens here, in one place (a copy on the real axis).
    static vec_t<T> AsT(const rvec_t& r)
    {
        vec_t<T> c(r.size());
        for (size_t g=0; g<r.size(); g++) c[g]=T(r[g]);
        return c;
    }
    //! The real-block tensor (\c DeltaPeriodicFaces' 3c-3 overload) and the transition pair's bodies: periodic-axis work,
    //! reached only through that base.
    const PairProjection<double>& RealTensor(const Orbital_DFT_IBS<double,dcmplx>& orb) const {return Tensor(itsO3R, itsPhiR, orb);}
    cvec_t        TransitionForwardImp(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                       const mat_t<dcmplx>& dD) const;
    mat_t<dcmplx> TransitionAdjointImp(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                       const cvec_t& v) const;
    //! ONE body for both \c Overlap3C overloads, each with its own typed cache.  R2.9(i) idiom: the
    //! accessor is const and the tensor is a lazily-built, geometry-fixed cache, so the maps are \c mutable.
    //! \brief ONE BLOCK'S \f$\Phi\f$, one of two realizations (V-ONE-OF-N; exactly one populated): the DENSE table, or the shared IMAGE
    //! record plus this block's phases, from which a contraction forms a point batch of rows when it needs one.  Read only by this class's
    //! own contraction bodies (they are the record's realization-aware clients); nothing above them sees which one is filled.
    template <class U> struct PhiRows
    {
        mat_t<U>                           dense;    //!< (npts x n), or empty
        std::shared_ptr<const PointImages> images;   //!< the shared record, or null
        cvec_t                             phases;   //!< this block's Bloch phase per image offset (with \c images)
        bool   IsDense()  const {return !images;}
        size_t NumRows()  const {return images ? images->NumPoints()    : dense.rows();}
        size_t NumCols()  const {return images ? images->NumFunctions() : dense.columns();}
    };
    template <class U> const PairProjection<U>& Tensor(const InstanceSlotMap<PairProjection<U>>& cache,
                                                   const InstanceSlotMap<PhiRows<U>>& phis,
                                                   const Orbital_1E_IBS<U>& orb) const;
    //! \f$\Phi_{gi}=\chi_i(r_g)\f$ -- HOW this class evaluates its own 3-centre integral, cached per block.
    //! Entirely private since 2026-08-23: a table of values is not an integral and has no business in an
    //! interface (user).  The three contractions below are what leaves.
    template <class U> const PhiRows<U>& Table(const InstanceSlotMap<PhiRows<U>>& cache, const Orbital_1E_IBS<U>& orb) const;
    //! Rows [lo,hi) of a block's \f$\Phi\f$ as a plain matrix, whichever realization holds it (a copy from a dense table; a phase sum
    //! from the images, bit-identical to the dense table's rows).
    template <class U> mat_t<U> Rows(const PhiRows<U>& P, size_t lo, size_t hi) const;
    //! The shared image record of \a src on my points, built once (threaded over point chunks) and cached per source.
    std::shared_ptr<const PointImages> Images(const BlochImageSource& src) const;
    //! \name THE GRADIENT TENSOR (GGA stage 2; \c DeltaRealFaces, real blocks only).
    //! \f$\partial_a\Phi_{gi}=\partial_a\chi_i(r_g)\f$, three dense (npts x n) tables per block, from the block's own \c VectorFunction::Gradient.
    //! THROWS on the images realization (stage 3) and on a block my factory evaluates itself (the atom's radial closure: it has no gradient
    //! yet, and the policy layer refuses a GGA on an atom).  The tensor's two closures are \c GradForward / \c GradAdjoint.
    //!@{
    using GradRows=std::array<mat_t<double>,3>;
    const PairProjection<double>& GradTensor(const Orbital_DFT_IBS<double,double>& orb) const;
    const GradRows& GradTable(const Orbital_1E_IBS<double>& orb) const;
    //! \f$[\nabla\rho]_{aN+g}=2\sum_i\partial_a\Phi_{gi}[\Phi D]_{gi}\f$ (D symmetric): one GEMM \f$\Phi D\f$ for all three components.
    rvec_t       GradForward(const mat_t<double>& P, const GradRows& G, const hmat_t<double>& D) const;
    //! \f$H=\Phi^TC+C^T\Phi\f$ with \f$C_{gi}=w_g\sum_a c_{aN+g}\partial_a\Phi_{gi}\f$: one GEMM.
    hmat_t<double> GradAdjoint(const mat_t<double>& P, const GradRows& G, const rvec_t& c) const;
    //!@}
    //! \f$\langle\chi_i|\sum_g c_g\delta_g|\chi_j\rangle=\Phi^\dagger\mathrm{diag}(w\,c)\Phi\f$
    template <class U> hmat_t<U> AdjointT(const PhiRows<U>& P, const rvec_t& c) const;
    template <class U> hmat_t<U> AdjointDenseT(const mat_t<U>& P, const rvec_t& c) const;
    template <class U> hmat_t<U> AdjointBatchedT(const PhiRows<U>& P, const rvec_t& c) const;
    //! \f$\langle\delta_g|\rho[D]\rangle/w_g=[\Phi D\Phi^\dagger]_{gg}\f$ -- the full quadratic form
    template <class U> rvec_t ForwardT(const PhiRows<U>& P, const hmat_t<U>& D) const;
    template <class U> rvec_t ForwardDenseT(const mat_t<U>& P, const hmat_t<U>& D) const;
    //! ...and the same thing for a caller holding \f$D=LL^\dagger\f$: \f$\sum_m|[\Phi L]_{gm}|^2\f$
    template <class U> rvec_t ForwardFactoredT(const PhiRows<U>& P, const mat_t<U>& L) const;
    // Per orbital BLOCK, keyed by the block's InstanceToken (S3c, user ruling 2026-11-14: the δ fit caches Φ itself).  Not by its
    // spatial sym_t: two blocks of equal symmetry (two basis sets, or one rebuilt at a new geometry against a surviving fit) are two
    // tables.  The spin channels share one block object, so they share its table.  A dead block's entries go on the next fill.
    // Thread-safe once-fills (the parallel Fock pass reaches them per block); a Tensor fill asks its Table, never the reverse.
    InstanceSlotMap<PhiRows<dcmplx>> itsPhi;    //!< Bloch blocks' Φ (npts x n), dense or images
    InstanceSlotMap<PhiRows<double>> itsPhiR;   //!< real TRIM blocks' Φ (3c-3)
    //! The shared image records, one per lattice basis behind the blocks (keyed by \c BlochImageSource::ImageSource), on the images route.
    InstanceSlotMap<std::shared_ptr<const PointImages>> itsImages;
    InstanceSlotMap<PairProjection<dcmplx>> itsO3;   //!< the Bloch blocks' 3-centre tensors, built on itsPhi
    InstanceSlotMap<PairProjection<double>> itsO3R;  //!< ...and the real TRIM blocks', on itsPhiR
    InstanceSlotMap<GradRows>               itsGradPhiR; //!< real blocks' ∂_aΦ (GGA stage 2), beside itsPhiR and keyed the same
    InstanceSlotMap<PairProjection<double>> itsGO3R;     //!< ...and their gradient tensors (\c GradTensor), on itsPhiR + itsGradPhiR
    FitQuadrature itsQuad;   //!< the mesh + its orbit fold + the Shubnikov spin tags
    RealPointValues itsRealValues;   //!< the factory's evaluation of a real block at my points, or empty (\c RealPointValues)
    ImageSourceOf   itsImageSourceOf;   //!< the factory's image source of a block, or empty (\c ImageSourceOf)
    PhiRealization  itsRealization=PhiRealization::Dense;   //!< set once, at setup, by \c ChoosePhiRealization
};

} //namespace

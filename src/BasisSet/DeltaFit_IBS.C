// File: BasisSet/DeltaFit_IBS.C  The delta (identity) v_xc fit basis over a real-space mesh.
//
// The OTHER representation a periodic Vxc fit basis can take (PlaneWaveFit_IBS being the {G} one): the
// delta basis, whose n_pts functions are the mesh's own delta functions.  Everything it knows is the
// quadrature -- points, weights, and the orbit fold the run's imposed ops gave that mesh -- because that
// IS what a delta basis is: W_g[h] = w_g delta_gh, a family of weight vectors with nothing left over.
//
// It exists so that CreateVxcFitBasisSet can make BOTH decisions -- which representation AND which points
// -- and return ONE object (doc/OpenWork.md, "separation of concerns in the XC terms").  Before this the
// mesh came back through a SECOND factory (CreateXCQuadrature) as a bare struct, and the Hamiltonian had
// to re-decide which term class to build in order to know which of the two answers to use; the fit basis
// and the grid were separate returns of one question.
//
// It ANSWERS, it does not hand out (user ruling 2026-08-22).  The mesh and its orbit fold are private:
// what leaves this class is an integral over its own functions or a projected/symmetrized coefficient
// vector -- arrays indexed by FUNCTION, in its own order.  No consumer knows whether these functions sit
// on an atom-centred Becke build or a uniform cell grid, which is exactly the fit/grid separation the
// design item is about.  Since 2026-08-22 NOTHING leaves: the mesh getter is gone with the Quadrature face.
//
// AND SINCE 2026-08-23 NOTHING SAYS "POINT" EITHER (user).  A delta basis is a family of FUNCTIONS, so it
// presents exactly the interface a Gaussian auxiliary basis presents: GetNumFunctions() for the count,
// Overlap(f) = <delta_g|f> for the projection, Integrals() = <delta_g|1> for the per-function integrals.
// The old NumPoints/Sample/Integrate trio described the QUADRATURE behind them, and read as correct only
// because here n_functions == n_points, so the wrong accessor returned the right number.
module;
#include <cassert>
#include <map>
#include <ostream>
#include <string>
#include <utility>
export module qchem.BasisSet.DeltaFit_IBS;
export import qchem.BasisSet.Fit_IBS;                    // cFIT_SF_Delta (the face) + FitQuadrature (the bundle)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.Symmetry;                                   // sym_t (the Bloch irrep)
import qchem.Types;                                      // dcmplx, rvec3_t, vec_t
import qchem.Mesh.Quadrature;                            // qcMesh::SiteIntegrals (the atomic-partition observable)
import qchem.Symmetry.Lattice_3D.Fold;                   // SymmetrizeValues / SymmetrizeValuesSigned (my orbit fold)
import qchem.Reporting;                                  // EmitFold -- this basis announces its own star-average

export namespace qchem::BasisSet
{

//! \brief The \f$\delta\f$-basis realization of the periodic \f$v_{xc}\f$ fit face: it holds the finished
//! XC quadrature and answers nothing else.
//!
//! Symmetry is the \f$\Gamma\f$ Bloch irrep (user ruling): the mesh is cell-periodic and carries no
//! crystal momentum of its own, exactly as \c PlaneWaveFit_IBS's fit block is built at \f$k=0\f$.
//! \c BasisSetID folds in the point count and the fold's orbit count -- the only things that distinguish
//! two of these -- but nothing keys an integral cache on a \f$\delta\f$ basis today: it has no 3-centre
//! tensor, because the quadrature route never forms one.
class DeltaFit_IBS
    : public virtual cFIT_SF_Delta            // FIT_SF_Delta<dcmplx> : cFIT_SF_ABS + Quadrature
    , public         IrrepBasisSetImp<dcmplx> // GetSymmetry/GetSymt/GetIrrep
{
public:
    //! \a q is the FINISHED quadrature (mesh + fold + Shubnikov tags) from the owning basis's
    //! \c CreateXCQuadrature, and \a sym the \f$\Gamma\f$ Bloch irrep.
    DeltaFit_IBS(FitQuadrature q, const sym_t& sym)
        : IrrepBasisSetImp<dcmplx>(sym)
        , itsQuad(std::move(q))
    {
        assert(itsQuad.mesh && itsQuad.mesh->size()>0 && "DeltaFit_IBS: a delta basis IS its mesh -- it cannot be empty");
        // The invariants of the bundle, checked HERE, once, where it becomes an object -- they used to be
        // re-checked by the consumer that received the loose struct.
        assert(itsQuad.fold.owner.empty()   || itsQuad.fold.owner.size()==itsQuad.mesh->size());
        assert(itsQuad.sigmas.empty()       || !itsQuad.fold.owner.empty());
        assert(itsQuad.flipFixed.empty()    || itsQuad.flipFixed.size()==itsQuad.mesh->size());
        // Providers self-report (user ruling 2026-08-16): the star-average this basis will apply to rho on
        // EVERY iteration announces itself at birth, armed or not.  nOps is known only for a magnetic fold
        // (sigmas run parallel to the ops); 0 = "not reported", NOT "not armed".
        qchem::report::EmitFold("XC mesh", itsQuad.sigmas.size(), itsQuad.mesh->size(),
                                itsQuad.fold.owner.empty() ? itsQuad.mesh->size() : itsQuad.fold.Reps(),
                                itsQuad.sigmas.empty() ? std::string() : std::string("magnetic (Shubnikov)"));
    }

    //! \copydoc BasisSet::FIT_SF_Delta::Overlap3C  (bodies in Imp/; ONE templated body serves both)
    const Projector3<double>& Overlap3C(const Orbital_1E_IBS<double>& orb) const override;
    const Projector3<dcmplx>& Overlap3C(const Orbital_1E_IBS<dcmplx>& orb) const override;

    // ---- the integrals over my own functions, one entry per FUNCTION (FIT_SF_ABS) ---------------------
    //! \copydoc BasisSet::FIT_SF_ABS::OverlapDiagonal
    //! \f$\langle\delta_g|\delta_g\rangle=w_g\f$: orthogonal, NOT orthonormal -- so a general fit
    //! through this basis divides by these, giving \f$c_g=w_g f_g/w_g=f(r_g)\f$, the point values.
    vec_t<dcmplx> OverlapDiagonal() const override {return AsComplex(itsQuad.mesh->Weights());}
    //! \copydoc BasisSet::FIT_SF_ABS::Charge
    //! \f$\langle\delta_g|1\rangle=\int\delta_g\,d^3r=w_g\f$.  Same numbers as \c OverlapDiagonal here
    //! and NOT the same question: that one is \f$\langle\delta_g|\delta_g\rangle\f$ (a metric), this is
    //! \f$\langle\delta_g|1\rangle\f$ (an integral).  They coincide because \f$\delta\f$ is idempotent
    //! under this quadrature -- a fact about this representation, not a shared definition.
    vec_t<dcmplx> Charge() const override {return AsComplex(itsQuad.mesh->Weights());}
    //! \copydoc BasisSet::FIT_SF_ABS::Overlap
    //! \f$\langle\delta_g|f\rangle=w_g f(r_g)\f$ -- I evaluate the field at MY points (the field's own
    //! bulk fast path does the work) and weight it.  The caller sees only a vector indexed by FUNCTION.
    vec_t<dcmplx> Overlap(const ScalarFunction<double>& f) const override
    {
        const rvec_t  v=f(itsQuad.mesh->Points());
        const rvec_t& w=itsQuad.mesh->Weights();
        assert(v.size()==w.size());
        vec_t<dcmplx> p(v.size());
        for (size_t g=0; g<v.size(); g++) p[g]=dcmplx(w[g]*v[g]);
        return p;
    }

    //! \copydoc BasisSet::FIT_SF_ABS::Symmetrize  (my mesh's orbit-mean projector)
    void Symmetrize(rvec_t& f) const override
    {
        if (itsQuad.fold.owner.empty()) return;              // free run: the projector is the identity
        Symmetry::Lattice_3D::SymmetrizeValues(itsQuad.fold, f);
    }
    //! \copydoc BasisSet::FIT_SF_ABS::SymmetrizeSpin
    //! The MAGNETIC case, and the only reason this override exists: with \f$\sigma\f$ tags the pair does
    //! NOT separate -- a Flip op maps \f$\rho_\uparrow\f$ onto \f$\rho_\downarrow\f$, not onto itself
    //! -- so \f$\rho\f$ takes the plain orbit mean while \f$m\f$ takes the \f$\chi\f$-signed one, with
    //! \f$m\f$ zeroed first at every function some Flip op fixes (where the exact projector annihilates
    //! it).  WITHOUT tags this falls through to the base's per-channel default, which is bit-identical to
    //! the branch that used to live here.
    void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const override
    {
        if (itsQuad.fold.owner.empty() || itsQuad.sigmas.empty())
            return FIT_SF_ABS<dcmplx>::SymmetrizeSpin(rho, m);   // free / grey: each channel on its own
        for (size_t g=0; g<itsQuad.flipFixed.size(); ++g) if (itsQuad.flipFixed[g]) m[g]=0.0;
        Symmetry::Lattice_3D::SymmetrizeValues      (itsQuad.fold, rho);
        Symmetry::Lattice_3D::SymmetrizeValuesSigned(itsQuad.fold, itsQuad.sigmas, m);
    }

    //! One \f$\delta\f$ function per mesh point.  ~100k of them: see the sizing warning on \c FIT_SF_Delta.
    size_t GetNumFunctions() const override {return itsQuad.mesh->size();}

    // NO op(r)/Gradient, and no stub either: since 2026-08-22 IrrepBasisSet does not promise pointwise
    // evaluation (doc/CleanupCandidates.md R1.0 step 1), so a basis whose functions are DISTRIBUTIONS
    // simply does not offer it.  This class used to carry two throwing overrides for exactly that reason
    // -- the interface demanded an answer it could not honestly give.  Nothing asks: the delta route's
    // fit coefficients ARE the field's values at the mesh points, and it has no expansion to plot, so it
    // does not implement FieldEvaluator either.

    std::string   Name      () const override {return "DeltaFit";}
    std::string   BasisSetID() const override
        {return Name()+"|npts="+std::to_string(itsQuad.mesh->size())
                      +"|orbits="+std::to_string(itsQuad.fold.owner.empty() ? itsQuad.mesh->size() : itsQuad.fold.Reps());}
    std::ostream& Write     (std::ostream& os) const override
        {return os << Name() << " fit IBS: " << itsQuad.mesh->size() << " delta functions on the XC mesh";}

private:
    //! The face is templated on the scalar of the ORBITAL blocks I contract against (Bloch => dcmplx), not
    //! on what my own functions are made of -- and a \f$\delta\f$ function and its weight are real.  So
    //! every integral over my own functions widens here, in one place.
    static vec_t<dcmplx> AsComplex(const rvec_t& r)
    {
        vec_t<dcmplx> c(r.size());
        for (size_t g=0; g<r.size(); g++) c[g]=dcmplx(r[g]);
        return c;
    }
    //! ONE body for both \c Overlap3C overloads, each with its own typed cache.  R2.9(i) idiom: the
    //! accessor is const and the tensor is a lazily-built, geometry-fixed cache, so the maps are \c mutable.
    template <class U> const Projector3<U>& Tensor(std::map<Irrep,Projector3<U>>& cache,
                                                   std::map<Irrep,mat_t<U>>& phis,
                                                   const Orbital_1E_IBS<U>& orb) const;
    //! \f$\Phi_{gi}=\chi_i(r_g)\f$ -- HOW this class evaluates its own 3-centre integral, cached per block.
    //! Entirely private since 2026-08-23: a table of values is not an integral and has no business in an
    //! interface (user).  The three contractions below are what leaves.
    template <class U> const mat_t<U>& Table(std::map<Irrep,mat_t<U>>& cache, const Orbital_1E_IBS<U>& orb) const;
    //! \f$\langle\chi_i|\sum_g c_g\delta_g|\chi_j\rangle=\Phi^\dagger\mathrm{diag}(w\,c)\Phi\f$
    template <class U> hmat_t<U> AdjointT(const mat_t<U>& P, const rvec_t& c) const;
    //! \f$\langle\delta_g|\rho[D]\rangle/w_g=[\Phi D\Phi^\dagger]_{gg}\f$ -- the full quadratic form
    template <class U> rvec_t ForwardT(const mat_t<U>& P, const hmat_t<U>& D) const;
    //! ...and the same thing for a caller holding \f$D=LL^\dagger\f$: \f$\sum_m|[\Phi L]_{gm}|^2\f$
    template <class U> rvec_t ForwardFactoredT(const mat_t<U>& P, const mat_t<U>& L) const;
    mutable std::map<Irrep,mat_t<dcmplx>> itsPhi;    //!< Bloch blocks' tables (npts x n)
    mutable std::map<Irrep,mat_t<double>> itsPhiR;   //!< real TRIM blocks' tables (disjoint irreps -- 3c-3)
    mutable std::map<Irrep,Projector3<dcmplx>> itsO3;   //!< the Bloch blocks' 3-centre tensors
    mutable std::map<Irrep,Projector3<double>> itsO3R;  //!< ...and the real TRIM blocks'
    FitQuadrature itsQuad;   //!< the mesh + its orbit fold + the Shubnikov spin tags
};

} //namespace

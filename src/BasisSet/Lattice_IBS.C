// File: BasisSet/Lattice_IBS.C  The LATTICE spec tier: the G=T concepts + the evaluator-templated IBS mixins.
//
// doc/BasisSetTaxonomyPlan.md §1.7 -- the four tiers.  This module is the SPEC tier for the translation
// group T: it says what an engine must answer to drive one Bloch block (the concepts) and turns a
// conforming engine into the role faces (the mixins).  It is tagged by the GROUP, not by a family: a
// plane-wave engine and a Gaussian engine (GPW) both satisfy it, which is why it lives HERE, in the core,
// BELOW every engine library -- it imports only the core faces.  (Before 1a0 it sat in qcLattice_BS as
// qchem.BasisSet.Lattice.IBS with the concepts inside the PW engine module; a spec that names a family
// is a defect, plan §1.7.)
//
// The lattice analog of BasisSet/Radial/IrrepBasisSet.C (Radial::Orbital_1E_IBS<E> etc.): the interface
// virtuals whose bodies are pure grid geometry are basis-agnostic and driven by the evaluator, so they
// live here once, templated on the evaluator E (constrained by isLattice_1E_Evaluator /
// isLattice_DFT_Evaluator), and a concrete lattice basis (PlaneWave_IBS, its auxiliary density-fit basis,
// GPW_IBS, APW/LAPW later) reuses them by instantiating the mixins with its own evaluator.
//
// As on the molecule/atom side, each mixin `dynamic_cast<const E&>(*this)`s (Cast()) to reach the
// evaluator base subobject of the final IBS (the IBS IS-A E -- a sibling base, so cross-cast RTTI, which
// is why every lattice evaluator has a polymorphic dtor).
//
// The concepts are purely STRUCTURAL (no derived_from): the spec knows no engine base class.  Check them
// where engine meets spec -- a static_assert beside each concrete IBS -- never inside an engine.
module;
#include <cassert>
#include <complex>      // std::real/std::imag (the ToScalar exact narrow)
#include <concepts>     // std::same_as (the concepts)
#include <functional>
#include <type_traits>  // std::is_same_v (ToScalar's identity branch)
export module qchem.BasisSet.Lattice_IBS;
import qchem.VectorFunction;                          // VectorFunction<T> (the pointwise-loop default of the point-SET op)
import qchem.BasisSet.IrrepBasisSet;                  // IrrepBasisSet<T> (op()(r), GetNumFunctions)
import qchem.BasisSet.Orbital_1E_IBS;                 // Orbital_1E_IBS<T> (MakeOverlap/MakeKinetic/MakeNuclear)
import qchem.BasisSet.Orbital_DFT_IBS;                // Orbital_DFT_IBS<T,dcmplx> (MakeRepulsion3C/MakeOverlap3C) + PairProjection<dcmplx>
import qchem.Structure;                               // Structure (MakeNuclear arg)
import qchem.Types;                                   // cvec_t, cvec3vec_t, chmat_t, rvec3_t

export namespace qchem::BasisSet::Lattice
{

//====================================================================================================
//  THE SPEC: what an engine must answer to drive a Bloch block.  Named for the G it serves (T), never
//  for a family -- a Gaussian engine (GPW) satisfies it without inheriting from the plane-wave one.
//====================================================================================================

//! \brief The lattice 1E evaluator concept the \c Orbital_1E_IBS mixin templates against (mirrors the
//! atom \c is1E_Evaluator).  Bloch-block evaluation (size/Eval/EvalGradient, complex: the block carries a
//! k-label) + the one-electron matrices (overlap/kinetic/nuclear) as \c chmat_t.  Models today:
//! \c PW_Evaluator (plane waves) and \c GPW_Evaluator (Bloch sums of Gaussians).
template <class E> concept isLattice_1E_Evaluator = requires (const E e, const rvec3_t& r, const Structure* cl)
{
    {e.size()             } -> std::same_as<size_t>;
    {e.Eval(r)            } -> std::same_as<cvec_t>;
    {e.EvalGradient(r)    } -> std::same_as<cvec3vec_t>;
    {e.OverlapMatrix()    } -> std::same_as<chmat_t>;
    {e.KineticMatrix()    } -> std::same_as<chmat_t>;
    {e.NuclearMatrix(cl)  } -> std::same_as<chmat_t>;
};

//! \brief The lattice DFT evaluator concept the \c Orbital_DFT_IBS mixin templates against (mirrors the
//! atom \c isDFT_Evaluator): the D-free reciprocal-space 3-centre tensors on top of the 1E tier.
//! (The potential->matrix bridge by G-vector LOOKUP that the old \c isPW_DFT_Evaluator also demanded is
//! NOT here: it is the PW fit-family pairing leaking into the G spec -- plan §1.7 -- and no mixin ever
//! consumed it.  It belongs to the fit side of \c Orbital_DFT_IBS<T,TFit>.)
template <class E> concept isLattice_DFT_Evaluator = isLattice_1E_Evaluator<E> && requires (const E e)
{
    {e.Repulsion3CTensor()} -> std::same_as<PairProjection<dcmplx>>;
    {e.Overlap3CTensor()  } -> std::same_as<PairProjection<dcmplx>>;
};

//====================================================================================================
//  EXACT NARROWING (doc/RealComplexPlan.md Step 3).  A TRIM block's Bloch matrices/vectors are real
//  BITWISE, not merely small-imaginary: Step 0 made the phases exactly ±1 (BlochPhase's parity form),
//  so every imaginary part is exactly 0.0 and the narrow is an ASSERTED fact, never a tolerance
//  (gate: GPW.TRIM_BlochMatricesAreExactlyReal).  ToScalar<T> is the identity for T=dcmplx, so the
//  mixin bodies below stay SINGLE-SOURCE across both instantiations.
//====================================================================================================
template <class T> T ToScalar(const dcmplx& z)
{
    if constexpr (std::is_same_v<T,dcmplx>) return z;
    else { assert(std::imag(z)==0.0 && "TRIM narrow: imaginary part must be EXACTLY zero (Step 0)"); return std::real(z); }
}
template <class T> vec_t<T> ToScalar(const cvec_t& v)
{
    if constexpr (std::is_same_v<T,dcmplx>) return v;
    else { vec_t<T> r(v.size()); for (size_t i=0;i<v.size();i++) r[i]=ToScalar<T>(v[i]); return r; }
}
template <class T> vec3vec_t<T> ToScalar(const cvec3vec_t& v)
{
    if constexpr (std::is_same_v<T,dcmplx>) return v;
    else
    {
        vec3vec_t<T> r(v.size());
        for (size_t i=0;i<v.size();i++) r[i]=vec3_t<T>(ToScalar<T>(v[i].x),ToScalar<T>(v[i].y),ToScalar<T>(v[i].z));
        return r;
    }
}
template <class T> mat_t<T> ToScalar(const mat_t<dcmplx>& m)
{
    if constexpr (std::is_same_v<T,dcmplx>) return m;
    else
    {
        mat_t<T> r(m.rows(), m.columns());
        for (size_t i=0;i<m.rows();i++) for (size_t j=0;j<m.columns();j++) r(i,j)=ToScalar<T>(m(i,j));
        return r;
    }
}
template <class T> hmat_t<T> ToScalar(const chmat_t& m)
{
    if constexpr (std::is_same_v<T,dcmplx>) return m;
    else
    {
        hmat_t<T> r(m.rows());
        for (size_t i=0;i<m.rows();i++)
            for (size_t j=i;j<m.columns();j++) r(i,j)=ToScalar<T>(m(i,j));
        return r;
    }
}

//====================================================================================================
//  THE MIXINS: a conforming engine -> the role faces.  Same names as the atom tier (Radial::Orbital_1E_IBS<E>),
//  qualified BasisSet::X where the core face of the same name is meant.
//====================================================================================================

// --- Shared tier: the IrrepBasisSet<T> evaluation + sizing that BOTH the orbital and the auxiliary
// (density-fit) lattice basis reuse, forwarded to the evaluator.  A cFIT_CD_ABS needs nothing more.
// T = the ORBITAL scalar (doc/RealComplexPlan.md Step 3), defaulted to dcmplx so every pre-Step-3 use
// (PlaneWave_IBS, the fit bases) is unchanged; a real TRIM block instantiates <E,double> and the
// evaluator's exactly-real results are ASSERT-narrowed by ToScalar (bodies single-source).
template <class E, class T=dcmplx> requires isLattice_1E_Evaluator<E>
class Irrep_IBS
    : public virtual BasisSet::Evaluatable_IBS<T>   // this mixin's whole job is op(r), so it makes the promise
{
public:
    virtual size_t       GetNumFunctions()          const override {return Cast().size();}
    virtual vec_t<T>     operator()(const rvec3_t& r) const override {return ToScalar<T>(Cast().Eval(r));}
    // d/dR, forward (GeometryOptimizationPlan §13.2): the gradient source awaits the user's ruling, so every family refuses for now.
    virtual bool CanGiveJets() const override {return false;}
    virtual Outcome<size_t,AD::NoGradient> JetsOn(const rvec3vec_t&, size_t, const Structure*, BasisSet::PointJets<T>&) const override
        {return BasisSet::NoGeometryDerivativeYet<size_t>("a lattice block", "function jets");}
    //! Batch: route to the evaluator's point-SET Bloch sum where it has one (GPW -- it pushes the image
    //! sum down into the molecular seam, so a transformed basis transforms once per POINT, not per
    //! image).  An evaluator without one (PW) keeps VectorFunction's pointwise default; the concept is
    //! deliberately NOT widened to demand it.
    virtual mat_t<T>     operator()(const rvec3vec_t& rs) const override
    {
        if constexpr (requires (const E& e) { e.EvalMany(rs); })
            return ToScalar<T>(Cast().EvalMany(rs));
        else
            return this->VectorFunction<T>::operator()(rs);
    }
    virtual vec3vec_t<T> Gradient  (const rvec3_t& r) const override {return ToScalar<T>(Cast().EvalGradient(r));}
protected:
    const E& Cast() const {return dynamic_cast<const E&>(*this);}
};

// --- Orbital 1E tier: the one-electron matrices (overlap/kinetic/nuclear), on top of the shared tier.
// Used by the ORBITAL bases (PlaneWave_IBS, GPW_IBS); an auxiliary fit basis does not carry these.
template <class E, class T=dcmplx> requires isLattice_1E_Evaluator<E>
class Orbital_1E_IBS
    : public Irrep_IBS<E,T>
    , public virtual BasisSet::Orbital_1E_IBS<T>
{
    using Irrep_IBS<E,T>::Cast;
public:
    virtual hmat_t<T> MakeOverlap() const override {return ToScalar<T>(Cast().OverlapMatrix());}
    virtual hmat_t<T> MakeKinetic() const override {return ToScalar<T>(Cast().KineticMatrix());}
    virtual hmat_t<T> MakeNuclear(const Structure* cl) const override {return ToScalar<T>(Cast().NuclearMatrix(cl));}
    // d/dR and d/dε (GeometryOptimizationPlan §13.2; doc/Records/CellHalf_Design.md §5).  An engine that answers the lattice-sum one-electron
    // pullbacks (GPW's Bloch Gaussians: libcint per image, over the energy series' own image list) gives the overlap and kinetic gradients AND
    // the image virials, one walk each; any other engine refuses.  PlaneWave_IBS overrides the pullbacks with exact zeros (its functions do
    // not move with the atoms).  A real TRIM block (T = double) hands its adjoint over widened: the Bloch sums at a TRIM k are real.
    virtual bool CanPullbackOverlap() const override {if constexpr (kOneElectronPullback) return Cast().CanOneElectronPullback(); else return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<T>& W, const Structure* st) const override
        {return PullbackOf(0, W, st, "overlap pullback");}
    virtual bool CanPullbackKinetic() const override {if constexpr (kOneElectronPullback) return Cast().CanOneElectronPullback(); else return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeKineticPullback(const hmat_t<T>& D, const Structure* st) const override
        {return PullbackOf(1, D, st, "kinetic pullback");}
    virtual bool CanVirialOverlap() const override {if constexpr (kOneElectronPullback) return Cast().CanOneElectronPullback(); else return false;}
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeOverlapVirial(const hmat_t<T>& W, const Structure* st) const override
        {return VirialOf(0, W, st, "overlap virial");}
    virtual bool CanVirialKinetic() const override {if constexpr (kOneElectronPullback) return Cast().CanOneElectronPullback(); else return false;}
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeKineticVirial(const hmat_t<T>& D, const Structure* st) const override
        {return VirialOf(1, D, st, "kinetic virial");}
    virtual bool CanPullbackNuclear() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeNuclearPullback(const hmat_t<T>&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a lattice block", "nuclear pullback");}
private:
    static constexpr bool kOneElectronPullback = requires (const E& e, const chmat_t& A, const Structure* st)
        { e.CanOneElectronPullback(); e.OneElectronPullback(0, A, st); };
    static chmat_t Widen(const hmat_t<T>& A)
    {
        if constexpr (std::is_same_v<T,dcmplx>) return A;
        else { chmat_t c(A.rows()); for (size_t i=0;i<A.rows();i++) for (size_t j=i;j<A.rows();j++) c(i,j)=dcmplx(A(i,j),0.0); return c; }
    }
    Outcome<rvec3vec_t,AD::NoGradient> PullbackOf(int which, const hmat_t<T>& A, const Structure* st, const char* what) const
    {
        if constexpr (kOneElectronPullback)
        {
            const auto r=Cast().OneElectronPullback(which, Widen(A), st);
            if (!r) return Outcome<rvec3vec_t,AD::NoGradient>::Fail(r.Error());
            return Outcome<rvec3vec_t,AD::NoGradient>::Ok(r->gradient);
        }
        else return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a lattice block", what);
    }
    Outcome<rmat3d_t,AD::NoGradient> VirialOf(int which, const hmat_t<T>& A, const Structure* st, const char* what) const
    {
        if constexpr (kOneElectronPullback)
        {
            const auto r=Cast().OneElectronPullback(which, Widen(A), st);
            if (!r) return Outcome<rmat3d_t,AD::NoGradient>::Fail(r.Error());
            return Outcome<rmat3d_t,AD::NoGradient>::Ok(r->virial);
        }
        else return BasisSet::NoGeometryDerivativeYet<rmat3d_t>("a lattice block", what);
    }
};

// --- Orbital DFT tier: the D-free reciprocal-space 3-centre tensors, forwarded to the evaluator.  Supplies
// the Orbital_DFT_IBS<T,dcmplx> one-time builds (the cached Repulsion3C/Overlap3C accessors call these).  The
// fit-basis arg is the delta support's declared cover (orbital-{G} intrinsic today), so it is not threaded
// to the evaluator yet -- GPW, whose density fit-grid does matter, threads it in its own override.
// NB the tensors follow TFit==dcmplx for BOTH orbital scalars (V1.1(ii): complex fit functions make <ab|c>
// complex even with real orbitals), so this tier needs NO narrowing -- exactly the split the
// Orbital_DFT_IBS<T,TFit> two-axis face prepared.
template <class E, class T=dcmplx> requires isLattice_DFT_Evaluator<E>
class Orbital_DFT_IBS
    : public virtual BasisSet::Orbital_DFT_IBS<T,dcmplx>
{
protected:
    virtual PairProjection<dcmplx> MakeRepulsion3C(const cFIT_CD_ABS&) const override {return Cast().Repulsion3CTensor();}
    virtual PairProjection<dcmplx> MakeOverlap3C  (const cFIT_SF_ABS&) const override {return Cast().Overlap3CTensor();}
    const E& Cast() const {return dynamic_cast<const E&>(*this);}
public:
    // d/dR (GeometryOptimizationPlan §13.2): the reciprocal-space Hartree has no 3-centre pullback (its force is the PW structure-factor
    // row in OpenWork §2, not this face).
    virtual bool CanPullbackRepulsion3C(const cFIT_CD_ABS&) const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsion3CPullback(const cFIT_CD_ABS&, const hmat_t<T>&, const cvec_t&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a lattice block", "3-centre repulsion pullback");}
};

} //namespace

// File: BasisSet/Lattice_3D/IrrepBasisSet.C  Lattice-generic, evaluator-templated plane-wave IBS mixins.
//
// The plane-wave analog of BasisSet/Molecule/IrrepBasisSet.C (module qchem.BasisSet.Molecule.IBS): the
// interface virtuals whose bodies are pure grid geometry are basis-agnostic and driven by the evaluator,
// so they live here once, templated on the evaluator E (constrained by isPW_1E_Evaluator / isPW_DFT_Evaluator), and a concrete
// plane-wave basis (the orbital PlaneWave_IBS, its auxiliary density-fit basis, and APW/LAPW later) reuses
// them by instantiating the mixins with its own evaluator.
//
// As on the molecule/atom side, each mixin `dynamic_cast<const E&>(*this)`s (Cast()) to reach the
// evaluator base subobject of the final IBS (the IBS IS-A E -- a sibling base, so cross-cast RTTI, which
// is why PW_Evaluator has a polymorphic dtor).
module;
#include <cassert>
#include <complex>      // std::real/std::imag (the ToScalar exact narrow)
#include <functional>
#include <type_traits>  // std::is_same_v (ToScalar's identity branch)
export module qchem.BasisSet.Lattice_3D.IBS;
import qchem.VectorFunction;                          // VectorFunction<T> (the pointwise-loop default of the point-SET op)
import qchem.BasisSet.IrrepBasisSet;                  // IrrepBasisSet<T> (op()(r), GetNumFunctions)
import qchem.BasisSet.Orbital_1E_IBS;                 // Orbital_1E_IBS<T> (MakeOverlap/MakeKinetic/MakeNuclear)
import qchem.BasisSet.Orbital_DFT_IBS;                    // Orbital_DFT_IBS<T,dcmplx> (MakeRepulsion3C/MakeOverlap3C) + Projector3<dcmplx>
import qchem.BasisSet.Fit_IBS;                        // cFIT_CD_ABS / cFIT_SF_ABS (the 3-centre fit-basis args)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;       // PW_Evaluator + isPW_1E_Evaluator / isPW_DFT_Evaluator
import qchem.Structure;                               // Structure (MakeNuclear arg)
import qchem.Types;                                   // cvec_t, cvec3vec_t, chmat_t, rvec3_t

export namespace qchem::BasisSet::Lattice_3D
{

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

// --- Shared tier: the IrrepBasisSet<T> evaluation + sizing that BOTH the orbital and the auxiliary
// (density-fit) plane-wave basis reuse, forwarded to the evaluator.  A cFIT_CD_ABS needs nothing more.
// T = the ORBITAL scalar (doc/RealComplexPlan.md Step 3), defaulted to dcmplx so every pre-Step-3 use
// (PlaneWave_IBS, the fit bases) is unchanged; a real TRIM block instantiates <E,double> and the
// evaluator's exactly-real results are ASSERT-narrowed by ToScalar (bodies single-source).
template <class E, class T=dcmplx> requires isPW_1E_Evaluator<E>
class EPW_Irrep_IBS
    : public virtual IrrepBasisSet<T>
{
public:
    virtual size_t       GetNumFunctions()          const override {return Cast().size();}
    virtual vec_t<T>     operator()(const rvec3_t& r) const override {return ToScalar<T>(Cast().Eval(r));}
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
// Used by the ORBITAL PlaneWave_IBS; an auxiliary fit basis does not carry these.
template <class E, class T=dcmplx> requires isPW_1E_Evaluator<E>
class EPW_Orbital1E_IBS
    : public EPW_Irrep_IBS<E,T>
    , public virtual Orbital_1E_IBS<T>
{
    using EPW_Irrep_IBS<E,T>::Cast;
public:
    virtual hmat_t<T> MakeOverlap() const override {return ToScalar<T>(Cast().OverlapMatrix());}
    virtual hmat_t<T> MakeKinetic() const override {return ToScalar<T>(Cast().KineticMatrix());}
    virtual hmat_t<T> MakeNuclear(const Structure* cl) const override {return ToScalar<T>(Cast().NuclearMatrix(cl));}
};

// --- Orbital DFT tier: the D-free reciprocal-space 3-centre tensors, forwarded to the evaluator.  Supplies
// the Orbital_DFT_IBS<T,dcmplx> one-time builds (the cached Repulsion3C/Overlap3C accessors call these).  The
// fit-basis arg is the delta support's declared cover (orbital-{G} intrinsic today), so it is not threaded
// to the evaluator yet -- GPW, whose density fit-grid does matter, threads it in its own override.
// NB the tensors follow TFit==dcmplx for BOTH orbital scalars (V1.1(ii): complex fit functions make <ab|c>
// complex even with real orbitals), so this tier needs NO narrowing -- exactly the split the
// Orbital_DFT_IBS<T,TFit> two-axis face prepared.
template <class E, class T=dcmplx> requires isPW_DFT_Evaluator<E>
class EPW_Orbital_DFT_IBS
    : public virtual Orbital_DFT_IBS<T,dcmplx>
{
protected:
    virtual Projector3<dcmplx> MakeRepulsion3C(const cFIT_CD_ABS&) const override {return Cast().Repulsion3CTensor();}
    virtual Projector3<dcmplx> MakeOverlap3C  (const cFIT_SF_ABS&) const override {return Cast().Overlap3CTensor();}
    const E& Cast() const {return dynamic_cast<const E&>(*this);}
};

} //namespace

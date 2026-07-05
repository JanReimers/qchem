// File: BasisSet/Lattice_3D/IrrepBasisSet.C  Lattice-generic, evaluator-templated plane-wave IBS mixins.
//
// The plane-wave analog of BasisSet/Molecule/IrrepBasisSet.C (module qchem.BasisSet.Molecule.IBS): the
// interface virtuals whose bodies are pure grid geometry are basis-agnostic and driven by the evaluator,
// so they live here once, templated on the evaluator E (constrained by isPW_Evaluator), and a concrete
// plane-wave basis (the orbital PlaneWave_IBS, its auxiliary density-fit basis, and APW/LAPW later) reuses
// them by instantiating the mixins with its own evaluator.
//
// As on the molecule/atom side, each mixin `dynamic_cast<const E&>(*this)`s (Cast()) to reach the
// evaluator base subobject of the final IBS (the IBS IS-A E -- a sibling base, so cross-cast RTTI, which
// is why PW_Evaluator has a polymorphic dtor).
module;
export module qchem.BasisSet.Lattice_3D.IBS;
import qchem.BasisSet.IrrepBasisSet;                  // IrrepBasisSet<dcmplx> (op()(r), GetNumFunctions)
import qchem.BasisSet.Orbital_1E_IBS;                 // Orbital_1E_IBS<dcmplx> (MakeOverlap/MakeKinetic)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;       // PW_Evaluator + isPW_Evaluator
import qchem.Types;                                   // cvec_t, cvec3vec_t, chmat_t, rvec3_t

export namespace qchem::BasisSet::Lattice_3D
{

// --- Shared tier: the IrrepBasisSet<dcmplx> evaluation + sizing that BOTH the orbital and the auxiliary
// (density-fit) plane-wave basis reuse, forwarded to the evaluator.  A cFIT_CD_ABS needs nothing more.
template <class E> requires isPW_Evaluator<E>
class EPW_Irrep_IBS
    : public virtual IrrepBasisSet<dcmplx>
{
public:
    virtual size_t     GetNumFunctions()          const override {return Cast().size();}
    virtual cvec_t     operator()(const rvec3_t& r) const override {return Cast().Eval(r);}
    virtual cvec3vec_t Gradient  (const rvec3_t& r) const override {return Cast().EvalGradient(r);}
protected:
    const E& Cast() const {return dynamic_cast<const E&>(*this);}
};

// --- Orbital 1E tier: the <p^2>/overlap building blocks (matrix-delivery), on top of the shared tier.
// Used by the ORBITAL PlaneWave_IBS; an auxiliary fit basis does not carry these.
template <class E> requires isPW_Evaluator<E>
class EPW_Orbital1E_IBS
    : public EPW_Irrep_IBS<E>
    , public virtual Orbital_1E_IBS<dcmplx>
{
    using EPW_Irrep_IBS<E>::Cast;
public:
    virtual chmat_t MakeOverlap() const override {return Cast().OverlapMatrix();}
    virtual chmat_t MakeKinetic() const override {return Cast().KineticMatrix();}
};

} //namespace

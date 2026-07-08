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
#include <functional>
export module qchem.BasisSet.Lattice_3D.IBS;
import qchem.BasisSet.IrrepBasisSet;                  // IrrepBasisSet<dcmplx> (op()(r), GetNumFunctions)
import qchem.BasisSet.Orbital_1E_IBS;                 // Orbital_1E_IBS<dcmplx> (MakeOverlap/MakeKinetic/MakeNuclear)
import qchem.BasisSet.Band_FT_IBS;                    // Band_FT_IBS (MakeRepulsion3C/MakeOverlap3C) + G_ERI3
import qchem.BasisSet.Fit_IBS;                        // cFIT_CD_ABS / cFIT_SF_ABS (the 3-centre fit-basis args)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;       // PW_Evaluator + isPW_1E_Evaluator / isPW_DFT_Evaluator
import qchem.Structure;                               // Structure (MakeNuclear arg)
import qchem.Types;                                   // cvec_t, cvec3vec_t, chmat_t, rvec3_t

export namespace qchem::BasisSet::Lattice_3D
{

// --- Shared tier: the IrrepBasisSet<dcmplx> evaluation + sizing that BOTH the orbital and the auxiliary
// (density-fit) plane-wave basis reuse, forwarded to the evaluator.  A cFIT_CD_ABS needs nothing more.
template <class E> requires isPW_1E_Evaluator<E>
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

// --- Orbital 1E tier: the one-electron matrices (overlap/kinetic/nuclear), on top of the shared tier.
// Used by the ORBITAL PlaneWave_IBS; an auxiliary fit basis does not carry these.
template <class E> requires isPW_1E_Evaluator<E>
class EPW_Orbital1E_IBS
    : public EPW_Irrep_IBS<E>
    , public virtual Orbital_1E_IBS<dcmplx>
{
    using EPW_Irrep_IBS<E>::Cast;
public:
    virtual chmat_t MakeOverlap() const override {return Cast().OverlapMatrix();}
    virtual chmat_t MakeKinetic() const override {return Cast().KineticMatrix();}
    virtual chmat_t MakeNuclear(const Structure* cl) const override {return Cast().NuclearMatrix(cl);}
};

// --- Orbital DFT tier: the D-free reciprocal-space 3-centre tensors, forwarded to the evaluator.  Supplies
// the Band_FT_IBS one-time builds (the cached Repulsion3C/Overlap3C accessors call these).  The fit-basis
// arg is the delta support's declared cover (orbital-{G} intrinsic today), so it is not threaded to the
// evaluator yet -- GPW, whose density fit-grid does matter, will thread it here.
template <class E> requires isPW_DFT_Evaluator<E>
class EPW_Orbital_DFT_IBS
    : public virtual Band_FT_IBS
{
public:
    //! The orbital-specific potential->KS-matrix bridge (Band_FT_IBS face): plane waves do the evaluator's
    //! Fourier lookup <i|V|j>=Vtilde(m_i-m_j).  GPW's Gaussian evaluator will instead grid-integrate here.
    virtual chmat_t MakeOverlap(const std::function<dcmplx(const ivec3_t&)>& Vt) const override
        {return Cast().OverlapMatrix(Vt);}
protected:
    virtual G_ERI3 MakeRepulsion3C(const cFIT_CD_ABS&) const override {return Cast().Repulsion3CTensor();}
    virtual G_ERI3 MakeOverlap3C  (const cFIT_SF_ABS&) const override {return Cast().Overlap3CTensor();}
    const E& Cast() const {return dynamic_cast<const E&>(*this);}
};

} //namespace

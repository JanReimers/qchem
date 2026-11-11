// File: Hamiltonian/Internal/IonIon.C  Ion-ion (nuclear-nuclear) repulsion ENERGY term.
//
// The T-templated collapse of the former Vnn (double, molecular) and PW_IonIon (dcmplx, plane-wave)
// terms -- they were the SAME energy-only term differing only in scalar type: each adds NO matrix
// contribution (a constant, independent of the electronic state) and asks the Structure its own ion-ion
// energy (Structure::IonIonEnergy): a direct Coulomb pair sum for a finite site set, an Ewald lattice sum for
// a periodic UnitCell (V-STRUCTURE-ISP F3).  The ion charge of each atom is supplied by a Z->charge callback: identity
// (itsZ) for the all-electron baseline, Z->Zion for a pseudopotential -- so the SAME term serves molecule,
// atom AND crystal, all-electron AND pseudopotential.  Mirrors how the two PP terms already share one model.
//
// Fully inline (no separate Imp TU): a template term must have its definition visible where it is
// instantiated -- IonIon<double> in the molecular Hamiltonians, IonIon<dcmplx> in the plane-wave ones.
module;
#include <type_traits>   // is_same_v -- the conditional real-block base (BuildBlocks)
#include <cassert>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <utility>
export module qchem.Hamiltonian.Internal.IonIon;
import qchem.Hamiltonian.Internal.Term;   // tStatic_HT<T> / tStatic_HT_Imp<T> (+ Energy/Types via re-export)
import qchem.Structure;                    // Structure::IonIonEnergy()/GetNumAtoms()
import qchem.Blaze;                        // blazem::zeroH<T> (the zero matrix contribution)
import qchem.Hamiltonian.Internal.GeometryGradient;   // the term's force face (GeometryOptimizationPlan §13.2)

export namespace qchem::Hamiltonian
{

//! \brief Ion-ion (nuclear-nuclear) repulsion energy \f$E_{nn}=\tfrac12\sum_{a\neq b} Z_a Z_b/|R_a-R_b|\f$
//! (finite) or its Ewald lattice sum (periodic), for any matrix element type \a T.
//!
//! A constant (density-independent) ENERGY term: it contributes NO Hamiltonian matrix (the ion-ion energy
//! does not depend on the electronic state), only \c EnergyBreakdown::Enn.  The pair-sum-vs-Ewald choice is
//! the structure's own (\c Structure::IonIonEnergy); \a T selects only the (zero) matrix
//! type (\c rsmat_t / \c chmat_t).  \c IonIon<double> is the molecular term, \c IonIon<dcmplx> the plane-wave
//! one.
template <class T> class IonIon
    : public virtual tStatic_HT<T>
    , private        tStatic_HT_Imp<T>
    , public         StaticRealBlockBase<T,IonIon<T>>   // real-block capability on the dcmplx instantiation (Step 3c)
    , public virtual tGeometryGradient_HT<T>            // its force (GeometryOptimizationPlan §13.2), after the private base (the cross-cast trap)
    , public virtual tStrainGradient_HT<T>              // its virial (CellHalf_Design §5; user ruling §6.2, 2026-11-10)
{
public:
    //! \copydoc HT_BlockOwner::BuildBlocks
    //! ⚠ THE CONDITIONAL-BASE CASE.  On \c <dcmplx> this term carries BOTH caches (its own and the
    //! real-block mixin's) and must prepare both; on \c <double> \c StaticRealBlockBase is \c NoRealBlock,
    //! which has no slots and no such method -- so the second call is compiled out rather than guarded at
    //! run time.  The ambiguity that forces this override on the complex instantiation is the shared
    //! \c HT_BlockOwner base doing its job (see it).
    virtual void BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override
    {
        tStatic_HT_Imp<T>::BuildBlocks(bs, cd);
        if constexpr (std::is_same_v<T,dcmplx>) StaticRealBlockBase<T,IonIon<T>>::BuildBlocks(bs, cd);
    }
    //! A constant energy: its (zero) block is the same for every spin channel, so ONE cached slot.
    virtual Spin CacheSpin(const Spin&) const override {return Spin::None;}
    typedef std::shared_ptr<const Structure> st_t;

    //! All-electron: the ion charge IS the true nuclear charge \c itsZ (identity map).
    IonIon(const st_t& st) : IonIon(st, [](int Z){return double(Z);}) {}

    //! Pseudopotential: \a zionOf maps an atom's true species \c Z to its ION CORE charge (the PP valence),
    //! so the atom's \c itsZ stays the true species while the Ewald/pair sum uses \c Zion.
    IonIon(const st_t& st, std::function<double(int)> zionOf)
        : tStatic_HT_Imp<T>()
        , theStructure(st)
        , itsZionOf(std::move(zionOf))
    {
        assert(theStructure && theStructure->GetNumAtoms()>0);
        assert(itsZionOf && "IonIon: a Z->ion-charge map is required");
        // EVALUATED ONCE, HERE.  E_nn depends only on the geometry and the ion charges, both fixed at
        // construction -- it is the most literally constant term in the Hamiltonian.  It used to be
        // recomputed on every GetEnergy call, i.e. a full EWALD LATTICE SUM per SCF iteration for a
        // periodic cell, to return the same number every time.  Computing it here also settles the
        // pair-sum-vs-Ewald choice at construction, where the structure is already known -- R2.16.  The
        // structure answers with its own sum (V-STRUCTURE-ISP F3: direct pairs, or Ewald for a UnitCell).
        itsEnn = theStructure->IonIonEnergy(itsZionOf);
    }

    virtual void GetEnergy(EnergyBreakdown& te, const tDM_CD<T>*) const override
    {
        te.Add("Enn", itsEnn, EnergyRole::Constant, 0.0);   // merges by name: a second Enn contributor adds, never clobbers (R1.2)
    }

    // --- the force (GeometryOptimizationPlan §13.2): the structure's own, beside its energy (Structure::IonIonForces) -- the pair sum's
    //     closed form for a finite site set, the Ewald sum's for a periodic cell -- with the same ion charges as the energy (itsZionOf). ---
    virtual std::string MissingGradient(const tbs_t<T>*, const Structure* st) const override
    {
        return st ? "" : "IonIon: no structure";
    }
    virtual void AddForces(rvec3vec_t& F, const tDM_CD<T>*, const Structure* st) const override
    {
        if (!st) throw std::logic_error("IonIon::AddForces: no structure");
        const rvec3vec_t f=st->IonIonForces(itsZionOf);
        if (f.size()!=F.size()) throw std::logic_error("IonIon::AddForces: the force array and the structure differ in atom count");
        for (size_t A=0;A<F.size();A++) F[A]+=f[A];
    }
    // --- the virial (CellHalf_Design §5): the structure's own, beside its force (Structure::IonIonVirial) -- the Ewald sum strained term by
    //     term for a periodic cell -- with the same ion charges.  A finite structure answers the pair sum's (Stress takes a cell, so a
    //     molecule never reaches a stress through it). ---
    virtual std::string MissingVirial(const tbs_t<T>*, const Structure* st) const override
    {
        return st ? "" : "IonIon: no structure";
    }
    virtual void AddVirial(rmat3d_t& V, const tDM_CD<T>*, const Structure* st) const override
    {
        if (!st) throw std::logic_error("IonIon::AddVirial: no structure");
        const rmat3d_t v=st->IonIonVirial(itsZionOf);
        for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) V(a,b)+=v(a,b);
    }
    virtual std::ostream& Write(std::ostream& os) const override
    {
        os << "    Ion-Ion energy E_nn = " << itsEnn << " Ha over " << theStructure->GetNumAtoms()
           << " ions (the structure's own sum: direct pairs, or Ewald for a cell)." << std::endl;
        return os;
    }

private:
    // No matrix contribution: a constant energy term, so the block is zero.
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>* bs, const Spin& s) const override {return MakeBlock<T>(bs,s);}
    //! ONE body for both block scalars (V1.35): zero, the native block and on <dcmplx> the real TRIM block (Static_HT_RealBlock_Of).
    template <class U> hmat_t<U> MakeBlock(const tobs_t<U>* bs, const Spin&) const
    {
        return blazem::zeroH<U>(bs->GetNumFunctions());
    }
    friend class Static_HT_RealBlock_Of<IonIon<T>>;

    st_t                       theStructure;
    std::function<double(int)> itsZionOf;   //!< Z -> ion core charge (identity for the all-electron baseline)
    double                     itsEnn=0.0;  //!< E_nn, evaluated ONCE in the ctor (geometry + charges are fixed)
};

} //namespace

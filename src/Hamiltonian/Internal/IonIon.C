// File: Hamiltonian/Internal/IonIon.C  Ion-ion (nuclear-nuclear) repulsion ENERGY term.
//
// The T-templated collapse of the former Vnn (double, molecular) and PW_IonIon (dcmplx, plane-wave)
// terms -- they were the SAME energy-only term differing only in scalar type: each adds NO matrix
// contribution (a constant, independent of the electronic state) and delegates to Ewald::NuclearRepulsion,
// which picks a direct Coulomb pair sum for a finite Structure and an Ewald lattice sum for a periodic cell
// (via Structure::isFinite()).  The ion charge of each atom is supplied by a Z->charge callback: identity
// (itsZ) for the all-electron baseline, Z->Zion for a pseudopotential -- so the SAME term serves molecule,
// atom AND crystal, all-electron AND pseudopotential.  Mirrors how the two PP terms already share one model.
//
// Fully inline (no separate Imp TU): a template term must have its definition visible where it is
// instantiated -- IonIon<double> in the molecular Hamiltonians, IonIon<dcmplx> in the plane-wave ones.
module;
#include <cassert>
#include <functional>
#include <iostream>
#include <memory>
#include <utility>
export module qchem.Hamiltonian.Internal.IonIon;
import qchem.Hamiltonian.Internal.Term;   // tStatic_HT<T> / tStatic_HT_Imp<T> (+ Energy/Types via re-export)
import qchem.Structure;                    // Structure::isFinite()/GetNumAtoms()
import qchem.Ewald;                        // NuclearRepulsion (pair sum for finite, Ewald for periodic)
import qchem.Blaze;                        // blazem::zeroH<T> (the zero matrix contribution)

export namespace qchem::Hamiltonian
{

//! \brief Ion-ion (nuclear-nuclear) repulsion energy \f$E_{nn}=\tfrac12\sum_{a\neq b} Z_a Z_b/|R_a-R_b|\f$
//! (finite) or its Ewald lattice sum (periodic), for any matrix element type \a T.
//!
//! A constant (density-independent) ENERGY term: it contributes NO Hamiltonian matrix (the ion-ion energy
//! does not depend on the electronic state), only \c EnergyBreakdown::Enn.  The pair-sum-vs-Ewald choice is
//! made downstream by \c NuclearRepulsion via \c Structure::isFinite(); \a T selects only the (zero) matrix
//! type (\c rsmat_t / \c chmat_t).  \c IonIon<double> is the molecular term, \c IonIon<dcmplx> the plane-wave
//! one.
template <class T> class IonIon
    : public virtual tStatic_HT<T>
    , private        tStatic_HT_Imp<T>
{
public:
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
    }

    virtual void GetEnergy(EnergyBreakdown& te, const tDM_CD<T>*) const override
    {
        // Direct pair sum for a finite molecule, Ewald lattice sum for a periodic cell (chosen by
        // isFinite()).  Charges come from the Z->ion map: itsZ all-electron, Zion for a pseudopotential.
        te.Enn = NuclearRepulsion(*theStructure, itsZionOf);
    }

    virtual std::ostream& Write(std::ostream& os) const override
    {
        size_t Na=theStructure->GetNumAtoms();
        if (theStructure->isFinite())
            os << "    Ion-Ion potential ZiZj/|Ri-Rj| with " << Na*(Na-1) << " ion pairs." << std::endl;
        else
            os << "    Ion-Ion (Ewald lattice sum) over " << Na << " ions per cell." << std::endl;
        return os;
    }

private:
    // No matrix contribution: a constant energy term, so the block is zero.
    virtual hmat_t<T> CalculateMatrix(const tobs_t<T>* bs, const Spin&) const override
    {
        return blazem::zeroH<T>(bs->GetNumFunctions());
    }

    st_t                       theStructure;
    std::function<double(int)> itsZionOf;   //!< Z -> ion core charge (identity for the all-electron baseline)
};

} //namespace

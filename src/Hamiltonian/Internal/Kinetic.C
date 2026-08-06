// File: Hamiltonian/Internal/Kinetic.C  Non-relativistic kinetic ENERGY term.
//
// The T-templated collapse of the former Kinetic (double, molecular) and PW_Kinetic (dcmplx,
// plane-wave) terms -- they were the SAME term differing only in scalar type.  Both computed
// 1/2<p^2> and contracted it against the density; the only real difference was that the plane-wave
// one called the UNCACHED MakeKinetic() while the molecular one used the cached Kinetic() accessor.
// That difference rested on a stale claim ("the symmetric cache is bypassed by the complex path"):
// Integrals_Kinetic<dcmplx> is explicitly instantiated (BasisSet/Imp/Orbital_1E_IBS.C), and its own
// comment there says the periodic dcmplx bases -- plane waves and GPW -- use these cached accessors
// too.  So the cached accessor serves both, and the periodic path now gets the cache for free
// (<p^2> is geometry-fixed, and the PW/GPW BasisSetID carries k, Ecut and nG).
//
// Fully inline (no separate Imp TU): a template term must have its definition visible where it is
// instantiated -- Kinetic<double> in the molecular Hamiltonians, Kinetic<dcmplx> in the plane-wave
// ones.  Mirrors IonIon<T> (the first of these collapses).
module;
#include <iostream>
export module qchem.Hamiltonian.Internal.Kinetic;
import qchem.Hamiltonian.Internal.Term;   // tStatic_HT<T> / tStatic_HT_Imp<T> (+ Energy/Types via re-export)
// Needed for blaze's scalar*matrix below.  A module that DEFINES a template using blaze operators must
// import this ITSELF: ADL at instantiation consults the template's DEFINITION context, so it is not
// enough that the instantiating TU imports it.  (Non-template code in an Imp TU never hits this -- it
// is why the old non-template Kinetic needed nothing.)  No <blaze/Math.h> include required.
import qchem.Blaze;

export namespace qchem::Hamiltonian
{

//! \brief Non-relativistic kinetic energy \f$T=-\tfrac12\nabla^2=\tfrac12\langle p^2\rangle\f$, for any
//! matrix element type \a T.
//!
//! Static (density-independent).  The basis supplies the \f$\langle p^2\rangle\f$ block with NO factor
//! \f$\tfrac12\f$ (see \c Integrals_Kinetic); this term applies the \f$\tfrac12\f$ -- i.e. the energy
//! boundary is where the factor lives.  \c Kinetic<double> is the molecular/atomic term,
//! \c Kinetic<dcmplx> the plane-wave/GPW one (diagonal in \f$|k+G|^2\f$ for a pure PW basis).
template <class T> class Kinetic
    : public virtual tStatic_HT<T>
    , private        tStatic_HT_Imp<T>
{
public:
    Kinetic() : tStatic_HT_Imp<T>() {}

    virtual void GetEnergy(EnergyBreakdown& te, const tDM_CD<T>* cd) const override
    {
        te.Kinetic += cd->DM_Contract(this);   // <T> = integral rho (1/2 p^2)
    }

    virtual std::ostream& Write(std::ostream& os) const override
    {
        return os << "    Kinetic energy 1/2<p^2>." << std::endl;
    }

private:
    virtual hmat_t<T> CalculateMatrix(const tobs_t<T>* bs, const Spin&) const override
    {
        return 0.5*bs->Kinetic();   // T = 1/2 <p^2>; Kinetic() is the cached <p^2> accessor
    }
};

} //namespace

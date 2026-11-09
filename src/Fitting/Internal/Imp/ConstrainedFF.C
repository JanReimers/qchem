// File: ConstrainedFF.C  Density (Coulomb-metric) fitter with the Dunlap charge constraint.
module;
#include <mutex>
#include <iostream>
#include <cassert>
#include <vector>
#include <stdexcept>
module qchem.Fitting.Internal.FunctionFitterImp;
import qchem.Fitting.Types;
import qchem.BasisSet.PairProjection;   // DensePairProjectionIntegrator + PairProjection (the borrowed tensor)
import qchem.Blaze;

namespace qchem::Fitting
{

//---------------------------------------------------------------------
//
//  Construction zone.  The Coulomb-metric fit replaces every overlap integral with a repulsion integral;
//  the initial guess carries the correct total charge (coeff[0] = 1/<f_0|1>).
//
template <class T> ConstrainedFF<T>::
ConstrainedFF(fbs_t& fbs, const vec_t<T>& theg)
    : itsBasisSet(fbs)
    , itsFitCoeff(fbs->GetNumFunctions(),0.0)
    , g  (theg)
    , gS (blazem::trans(g)*fbs->InvRepulsion())
    , gSg(gS*g)
{
    itsFitCoeff[0]=1.0/itsBasisSet->Charge()[0]; // wild guess with the correct total charge
}

template <class T> const DensePairProjectionIntegrator<T>&
ConstrainedFF<T>::Integrator(const BasisSet::Orbital_DFT_IBS<T,T>& orb) const
{
    std::lock_guard<std::mutex> lock(itsIntLock.m);   // a parallel Fock pass's blocks share this fitter (see DeltaScalarFitter::Integrator)
    const sym_t& id=orb.GetSymt();
    auto it=itsInt.find(id);
    if (it!=itsInt.end()) return it->second;
    // The basis's cached, D-free Coulomb tensor <ab|c> (built once, keyed by BasisSetID) -- borrowed by the
    // integrator, which is a VIEW of it; the fit functions' own <f_a|1> size its energy quadrature.
    const PairProjection<T>& R3=orb.Repulsion3C(*itsBasisSet);
    return itsInt.emplace(id, DensePairProjectionIntegrator<T>(R3, itsBasisSet->Charge())).first->second;
}

template <class T> const qcMesh::MatrixForward<double>&
ConstrainedFF<T>::Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    return Integrator(orb);
}

template <class T> void ConstrainedFF<T>::Constrain(T N)
{
    T lambda = (blazem::trans(g)*itsFitCoeff - N) / gSg;
    itsFitCoeff -= lambda * blazem::trans(gS);                // enforce g.c = N exactly
}

// The density's Coulomb projections (a matrix density's forward, or a seed's ProjectField) and its charge.  Robust / variational
// density fitting with a linear (charge) constraint, after
//   B. I. Dunlap, J. W. D. Connolly & J. R. Sabin, J. Chem. Phys. 71(8), 3396 (1979):
// the unconstrained Coulomb-metric fit c0 = G^-1 p, then g.c = N enforced exactly by one Lagrange correction,
//   c = c0 - lambda G^-1 g,  lambda = (g.c0 - N)/(g.G^-1 g) = (g.c0 - N)/gSg,  G^-1 g = trans(gS).
// It minimizes the Coulomb self-energy of the residual at exact charge, so the fitted Vee is variational (second order in the fit error).
template <class T> void ConstrainedFF<T>::DoFit(const Projections& rho, double charge)
{
    itsFitCoeff = itsBasisSet->InvRepulsion() * rho.Vector();   // c0 = G^-1 <f|1/r12|rho>
    Constrain(charge);
}

template <class T> hmat_t<T> ConstrainedFF<T>::Repulsion(const BasisSet::Orbital_DFT_IBS<T,T>& orb) const
{
    hmat_t<T> J=Integrator(orb).Adjoint(itsFitCoeff);   // THE ADJOINT HALF of the forward the density contracted D into
    assert(!blazem::isnan(J));
    return J;
}

template <class T> Projections ConstrainedFF<T>::ProjectField(const ScalarFunction<double>& f) const
{
    const auto* sf=dynamic_cast<const BasisSet::FIT_SF_NonOrtho*>(itsBasisSet.get());   // the overlap metric-solve face
    if (!sf) throw std::logic_error("ConstrainedFF::ProjectField: the fit basis has no overlap metric (FIT_SF_NonOrtho face)");
    const rvec_t c0=sf->InvOverlap()*sf->Overlap(f);                       // S^-1 <f|rho>
    return Projections(rvec_t(itsBasisSet->Repulsion()*c0));          // G c0: the Coulomb projections of that fit
}

//---------------------------------------------------------------------------
//
//  Fit-derived quantities the clients query (the "what's your repulsion with this basis?" side).
//
template <class T> double ConstrainedFF<T>::SelfRepulsion() const
{
    // The CROSS metric on itself (the form the old FitGetRepulsion used), so the energy keeps its bits.
    return 0.5*(blazem::trans(itsFitCoeff) * itsBasisSet->Repulsion(*itsBasisSet) * itsFitCoeff);   // 1/2 <fit|1/r12|fit>
}

//---------------------------------------------------------------------------
//
//  The geometry derivatives at fixed c (GeometryOptimizationPlan §13.2).  The Dunlap fit is variational, so a force needs only these: the
//  c-derivative of its energy is λ g with g = ∫f_k fixed (a function's integral does not move with it), and g·c = N, so c's response
//  enters only as λ dN/dR, which is ZERO: N is the fixed electron count (user ruling 2026-11-08; Hamiltonian/Internal/Imp/Vee_FitDunlap_Forces.C).
//
template <class T> bool ConstrainedFF<T>::CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<T,T>& orb) const
{
    return orb.CanPullbackRepulsion3C(*itsBasisSet);
}

template <class T> Outcome<rvec3vec_t,AD::NoGradient> ConstrainedFF<T>::MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<T,T>& orb, const hmat_t<T>& D,
                                                                                          const Structure* st) const
{
    // Repulsion(orb) = Σ_k c_k (ij|k) exactly (the integrator's Adjoint has no screen), so its pullback is the tensor's against D ⊗ c.
    return orb.MakeRepulsion3CPullback(*itsBasisSet, D, itsFitCoeff, st);
}

template <class T> bool ConstrainedFF<T>::CanPullbackSelfRepulsion() const
{
    return itsBasisSet->CanPullbackRepulsion();
}

template <class T> Outcome<rvec3vec_t,AD::NoGradient> ConstrainedFF<T>::MakeSelfRepulsionPullback(const Structure* st) const
{
    // SelfRepulsion = ½ cᵀ G c, so the adjoint of G is ½ c cᵀ (upper triangle; the basis face reads it symmetric).
    const vec_t<T>& c=itsFitCoeff;
    rsmat_t A(c.size());
    for (size_t k=0;k<c.size();k++) for (size_t l=k;l<c.size();l++) A(k,l)=0.5*c[k]*c[l];
    return itsBasisSet->MakeRepulsionPullback(A, st);
}

template <class T> double ConstrainedFF<T>::Integral() const
{
    const vec_t<T> q=itsBasisSet->Charge();   // <f_a|1>, by value (one signature for both fit faces)
    return blazem::trans(itsFitCoeff) * q;
}

template class ConstrainedFF<double>;

} //namespace

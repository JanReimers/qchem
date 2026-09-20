// File: ConstrainedFF.C  Density (Coulomb-metric) fitter with the Dunlap charge constraint.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.Fitting.Internal.FunctionFitterImp;
import qchem.Fitting.Types;
import qchem.Streamable;
import qchem.Blaze;
import qchem.BasisSet.Projector3;   // DenseProjector3Integrator + Projector3 (the borrowed tensor)

namespace qchem::Fitting
{

//---------------------------------------------------------------------
//
//  Construction zone.  The Coulomb-metric fit replaces every overlap integral with a repulsion integral;
//  the initial guess carries the correct total charge (coeff[0] = 1/<f_0|1>).
//
template <class T> ConstrainedFF<T>::ConstrainedFF()
    : Base()
    , g  ( )
    , gS ( )
    , gSg(0)
{}

template <class T> ConstrainedFF<T>::
ConstrainedFF(fbs_t& fbs, const vec_t<T>& theg)
    : Base(fbs)
    , g  (theg)
    , gS (blazem::trans(g)*fbs->InvRepulsion())
    , gSg(gS*g)
{
    this->itsFitCoeff[0]=1.0/this->itsBasisSet->Charge()[0]; // wild guess with the correct total charge
}

//--------------------------------------------------------------------------
//
//  DoFit:  c0 = the projection's own-metric unconstrained fit + the Dunlap charge constraint.  The metric is
//  the projection's business (a STRATEGY dispatched by polymorphism): a density with a matrix returns the
//  Coulomb-metric solve J^-1 <rho|c> (the ProjectedDensity_AO default); a matrix-free seed returns its
//  overlap-metric fit S^-1<f|rho> directly.  The fitter owns only the charge constraint below.
//
template <class T> void ConstrainedFF<T>::DoFitUnconstrained(const ProjectedDensity_AO& ffc)
{
    // Which metric face the projection HAS decides the route (V1.16: a face you have or do not, never a
    // default that guesses).  Both are "what can it do" cross-casts between abstract faces.
    if (auto* cm=dynamic_cast<const CoulombMetric_ProjectedDensity*>(&ffc))
    {
        // A matrix-carrying density: it contracts its own D into the forward THIS fitter vends, per block,
        // and the metric solve is ours -- c0 = J^-1 <rho|c>.  (R1.0q, 2026-09-19: the solve moved here
        // from the projection, which had no business owning the fit basis's metric.)
        this->itsFitCoeff = this->itsBasisSet->InvRepulsion() * cm->GetRepulsion3C(*this);
        return;
    }
    auto* om=dynamic_cast<const OverlapMetric_ProjectedDensity*>(&ffc);
    assert(om && "ConstrainedFF: a ProjectedDensity_AO carries either the Coulomb or the overlap metric face");
    this->itsFitCoeff = om->GetUnconstrainedFit(this->itsBasisSet.get());   // a seed's own S^-1 <f|rho>
}

template <class T> const DenseProjector3Integrator<T>&
ConstrainedFF<T>::Integrator(const BasisSet::Orbital_DFT_IBS<T,T>& orb) const
{
    const sym_t& id=orb.GetSymt();
    auto it=itsInt.find(id);
    if (it!=itsInt.end()) return it->second;
    // The basis's cached, D-free Coulomb tensor <ab|c> (built once, keyed by BasisSetID) -- borrowed by the
    // integrator, which is a VIEW of it; the fit functions' own <f_a|1> size its energy quadrature.
    const Projector3<T>& R3=orb.Repulsion3C(*this->itsBasisSet);
    return itsInt.emplace(id, DenseProjector3Integrator<T>(R3, this->itsBasisSet->Charge())).first->second;
}

template <class T> const qcMesh::MatrixForward<double>&
ConstrainedFF<T>::Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    return Integrator(orb);
}

template <class T> void ConstrainedFF<T>::DoFit(const ProjectedDensity<T>& pd)
{
    // NON-orthonormal (Gaussian) density fit: recover the AO projection face -- a sanctioned abstract->abstract
    // cross-cast, the paired-fitter half of the neutral ProjectedDensity<T> seam (the {G}-map projection is the
    // orthonormal sibling, consumed by the reciprocal-space fitter instead).
    const auto* ffcp = dynamic_cast<const ProjectedDensity_AO*>(&pd);
    assert(ffcp && "ConstrainedFF (non-ortho Gaussian density fit) requires a ProjectedDensity_AO projection");
    const ProjectedDensity_AO& ffc = *ffcp;
    // Robust / variational density fitting with a linear (charge) constraint, after
    //   B. I. Dunlap, J. W. D. Connolly & J. R. Sabin, J. Chem. Phys. 71(8), 3396 (1979).
    // Do the unconstrained Coulomb-metric fit  c0 = J^-1 b  (J = Coulomb/repulsion metric, b = <f_a|rho>),
    // then enforce  g.c = N  exactly (g_a = integral f_a, N = total charge) by one Lagrange correction:
    //   c = c0 - lambda J^-1 g,  lambda = (g.c0 - N)/(g.J^-1 g) = (g.c0 - N)/gSg,  J^-1 g = trans(gS).
    // Minimizes the Coulomb self-energy of the residual subject to exact charge, so the fitted Vee is
    // variational (error second order in the fit error) rather than relying on a post-hoc rescale.
    DoFitUnconstrained(ffc);                                        // c0 -> itsFitCoeff (unconstrained)
    T N      = ffc.FitGetConstraint();
    T lambda = (blazem::trans(g)*this->itsFitCoeff - N) / gSg;
    this->itsFitCoeff -= lambda * blazem::trans(gS);                // enforce g.c = N exactly
}

//---------------------------------------------------------------------------
//
//  Fit-derived quantities the clients query (the "what's your repulsion with this basis?" side).
//
template <class T> hmat_t<T> ConstrainedFF<T>::Repulsion(const robs_t<T>* bs) const
{
    // robs_t is the 1E base; the 3-centre tier is the DFT one -- the cross-cast the periodic fitters also make.
    const auto& dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<T,T>&>(*bs);
    // THE ADJOINT HALF, off the same object whose forward the density contracted D into (R1.0q).
    hmat_t<T> J=Integrator(dftbs).Adjoint(this->itsFitCoeff);
    assert(!blazem::isnan(J));
    return J;
}

template <class T> double ConstrainedFF<T>::FitGetRepulsion(const ConstrainedFF<T>* ffi) const
{
    return
        blazem::trans(this->itsFitCoeff) * this->itsBasisSet->Repulsion(*ffi->itsBasisSet.get()) *
        ffi->itsFitCoeff;
}

template <class T> double ConstrainedFF<T>::FitGetSelfRepulsion() const
{
    return FitGetRepulsion(this);   // <fit|1/r12|fit>
}

template <class T> double ConstrainedFF<T>::Integral() const
{
    const vec_t<T> q=this->itsBasisSet->Charge();   // <f_a|1>, by value (one signature for both fit faces)
    return blazem::trans(this->itsFitCoeff) * q;
}

template <class T> std::ostream& ConstrainedFF<T>::Write(std::ostream& os) const
{
    Base::Write(os);
    os << g << gS;
    return os;
}

template class ConstrainedFF<double>;

} //namespace

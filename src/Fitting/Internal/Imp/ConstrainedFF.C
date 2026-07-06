// File: ConstrainedFF.C  Density (Coulomb-metric) fitter with the Dunlap charge constraint.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.Fitting.Internal.FunctionFitterImp;
import qchem.Fitting.Types;
import qchem.Streamable;
import qchem.Blaze;

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
    this->itsFitCoeff = ffc.GetUnconstrainedFit(this->itsBasisSet.get());
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
    auto dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<T>*>(bs); // robs_t is the 1E base; need the 3-centre one
    assert(dftbs && "ConstrainedFF::Repulsion: Gaussian fitting needs an Orbital_DFT_IBS (3-centre) basis");
    const ERI3<T>& R3=dftbs->Repulsion3C(*this->itsBasisSet);
    hmat_t<T> J=blazem::zeroH<T>(bs->GetNumFunctions());
    size_t i=0;
    for (auto c:this->itsFitCoeff) J+=c*R3[i++];
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
    return blazem::trans(this->itsFitCoeff) * this->itsBasisSet->Charge();
}

template <class T> std::ostream& ConstrainedFF<T>::Write(std::ostream& os) const
{
    Base::Write(os);
    os << g << gS;
    return os;
}

template class ConstrainedFF<double>;

} //namespace

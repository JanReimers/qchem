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
//  DoFit:  c0 = J^-1 <rho|f>  (unconstrained Coulomb-metric solve) + the Dunlap charge constraint.
//
template <class T> void ConstrainedFF<T>::DoFitUnconstrained(const ProjectedDensity_AO& ffc)
{
    auto Jinv=this->itsBasisSet->InvRepulsion();
    this->itsFitCoeff=Jinv * ffc.GetRepulsion3C(this->itsBasisSet.get());
}

template <class T> void ConstrainedFF<T>::DoFit(const ProjectedDensity_AO& ffc)
{
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
template <class T> hmat_t<T> ConstrainedFF<T>::Repulsion(const obs_t<T>* bs) const
{
    auto dftbs=dynamic_cast<const BasisSet::Orbital_DFT_IBS<T>*>(bs); // obs_t is the 1E base; need the 3-centre one
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

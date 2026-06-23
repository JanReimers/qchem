// File: ConstrainedFF.C  General constrained fit.
module;
#include <iostream>
module qchem.FittedFunctionImp;
import qchem.Fitting.Types;
import qchem.Blaze;

namespace qchem::Fitting
{

template <class T> ConstrainedFF<T>::ConstrainedFF()
    : FittedFunctionImp<T>()
    , g  ( )
    , gS ( )
    , gSg(0)
{}

template <class T> ConstrainedFF<T>::
ConstrainedFF(bs_t& fbs, const vec_t<T>& theg, mesh_t&  m)
    : FittedFunctionImp<T>(fbs,m)
    , g  (theg)
    , gS (blazem::trans(g)*fbs->InvRepulsion())
    , gSg(gS*g)
{	
}

template <class T> void ConstrainedFF<T>::DoFit(const ScalarFFClient& ffc)
{
    FittedFunctionImp<T>::DoFitInternal(ffc);   // scalar (potential) fit: overlap metric, unconstrained
}
template <class T> void ConstrainedFF<T>::DoFit(const DensityFFClient& ffc)
{
    // Robust / variational density fitting with a linear (charge) constraint, after
    //   B. I. Dunlap, J. W. D. Connolly & J. R. Sabin, J. Chem. Phys. 71(8), 3396 (1979).
    // Do the unconstrained Coulomb-metric fit  c0 = J^-1 b  (J = Coulomb/repulsion metric, b = <f_a|rho>),
    // then enforce  g.c = N  exactly (g_a = integral f_a, N = total charge) by one Lagrange correction:
    //   c = c0 - lambda J^-1 g,  lambda = (g.c0 - N)/(g.J^-1 g) = (g.c0 - N)/gSg,  J^-1 g = trans(gS).
    // Minimizes the Coulomb self-energy of the residual subject to exact charge, so the fitted Vee is
    // variational (error second order in the fit error) rather than relying on a post-hoc rescale.
    FittedFunctionImp<T>::DoFitInternal(ffc);                       // c0 -> itsFitCoeff (unconstrained)
    T N      = ffc.FitGetConstraint();
    T lambda = (blazem::trans(g)*this->itsFitCoeff - N) / gSg;
    this->itsFitCoeff -= lambda * blazem::trans(gS);                // enforce g.c = N exactly
}

template <class T> std::ostream& ConstrainedFF<T>::Write(std::ostream& os) const
{
    FittedFunctionImp<T>::Write(os);
    os << g << gS;
    return os;
}


template class ConstrainedFF<double>;

} //namespace
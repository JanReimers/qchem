// File: ConstrainedFF.C  General constrained fit.
module;
#include <iostream>
#include "blaze/Math.h"
module qchem.FittedFunctionImp;
import qchem.FittedFunctionClient;
import qchem.Fit_IBS;


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
    , gS (blaze::trans(g)*fbs->InvRepulsion())
    , gSg(gS*g)
{	
}

template <class T> void ConstrainedFF<T>::DoFit(const ScalarFFClient& ffc)
{
    FittedFunctionImp<T>::DoFitInternal(ffc);
}
template <class T> void ConstrainedFF<T>::DoFit(const DensityFFClient& ffc)
{
    FittedFunctionImp<T>::DoFitInternal(ffc,ffc.FitGetConstraint());
}

template <class T> std::ostream& ConstrainedFF<T>::Write(std::ostream& os) const
{
    FittedFunctionImp<T>::Write(os);
    os << g << gS;
    return os;
}


template class ConstrainedFF<double>;

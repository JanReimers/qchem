// File: ConstrainedFF.C  General constrained fit.
module;
#include <iostream>
#include <cassert>
#include <vector>
module qchem.FittedFunctionImp;
import qchem.FittedFunctionClient;
import qchem.Fit_IBS;
import oml;


template <class T> ConstrainedFF<T>::ConstrainedFF()
    : FittedFunctionImp<T>()
    , g  ( )
    , gS ( )
    , gSg(0)
{}

template <class T> ConstrainedFF<T>::
ConstrainedFF(bs_t& fbs, const Vec& theg, mesh_t&  m)
    : FittedFunctionImp<T>(fbs,m)
    , g  (theg)
    , gS (g*fbs->InvRepulsion(itsLAParams))
    , gSg(gS*g)
{	
}

template <class T> double ConstrainedFF<T>::DoFit(const ScalarFFClient& ffc)
{
    return FittedFunctionImp<T>::DoFitInternal(ffc);
}
template <class T> double ConstrainedFF<T>::DoFit(const DensityFFClient& ffc)
{
    return FittedFunctionImp<T>::DoFitInternal(ffc,ffc.FitGetConstraint());
}

template <class T> std::ostream& ConstrainedFF<T>::Write(std::ostream& os) const
{
    FittedFunctionImp<T>::Write(os);
    os << g << gS;
    return os;
}


template class ConstrainedFF<double>;

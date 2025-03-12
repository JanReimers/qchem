// File: ConstrainedFF.C  General constrained fit.

#include "Imp/Fitting/ConstrainedFF.H"
#include <FittedFunctionClient.H>
#include <BasisSet.H>
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>


template <class T> ConstrainedFF<T>::ConstrainedFF()
    : FittedFunctionImp<T>()
    , g  ( )
    , gS ( )
    , gSg(0)
{}

template <class T> ConstrainedFF<T>::
ConstrainedFF(bs_t& theFitBasisSet, const Vec& theg, mesh_t&  m)
    : FittedFunctionImp<T>(theFitBasisSet,m)
    , g  (theg)
    , gS (g*itsInvRepl)
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
    if (StreamableObject::Binary())
    {
        BinaryWrite(gSg,os);
    }
    else if(StreamableObject::Ascii())
    {
        os << gSg << " ";
    }
    return os;
}

template <class T> std::istream& ConstrainedFF<T>::Read (std::istream& is)
{
    FittedFunctionImp<T>::Read(is);
    is >> g >> gS;
    if (StreamableObject::Binary())
    {
        BinaryRead(gSg,is);
    }
    else
    {
        is >> gSg;
        assert(is.get() == ' ');
    }
    return is;
}

template class ConstrainedFF<double>;

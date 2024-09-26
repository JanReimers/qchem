// File: ConstrainedFF.C  General constrained fit.



#include "FunctionsImp/ConstrainedFF.H"
#include "Functions/FittedFunctionClient.H"
#include "BasisSet/TBasisSet.H"
#include "BasisSet/IntegralDataBase.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
//#include "oml/vector_io.h"
#include <iostream>
#include <cassert>


template <class T> ConstrainedFF<T>::ConstrainedFF()
    : FittedFunctionImplementation<T>()
    , g  ( )
    , gS ( )
    , gSg(0)
{}

template <class T> ConstrainedFF<T>::
ConstrainedFF(const rc_ptr<BasisSet>& theFitBasisSet, const Vec& theg, bool CDfit)
    : FittedFunctionImplementation<T>(theFitBasisSet,CDfit)
    , g  (theg)
    , gS (g*GetInverseOverlap())
    , gSg(gS*g)
{	
}

template <class T> double ConstrainedFF<T>::DoFit(const FittedFunctionClient& ffc)
{
    return FittedFunctionImplementation<T>::DoFit(ffc.FitGetConstraint(),ffc);
}

template <class T> double ConstrainedFF<T>::
DoFit(double constraint, const Vec& overlap)
{
    double lam= gSg>0 ? (constraint-gS*overlap)/gSg : 0;
    /*	cout << "ConstrainedFF<T>::DoFit" << lam << " " << gSg << " " << constraint << endl
    		<< gS << endl << overlap << endl << g << std::endl;*/
    FittedFunctionImplementation<T>::SetFitCoeff(GetInverseOverlap()*(overlap+lam*g));
    return lam;
}

template <class T> std::ostream& ConstrainedFF<T>::Write(std::ostream& os) const
{
    FittedFunctionImplementation<T>::Write(os);
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
    FittedFunctionImplementation<T>::Read(is);
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

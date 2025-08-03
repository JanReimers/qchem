// File: TBasisSetImplementation.C
module;
#include <tuple>
#include <iostream>
#include <vector>
import qchem.LAParams;
import qchem.LASolver;

module qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisFunction;
import oml;

LAParams DefaultLAP({qchem::Lapack,qchem::SVD,1e-10,1e-12});
//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> Orbital_IBS_Common1<T>::Orbital_IBS_Common1()
    : itsLAParams      (DefaultLAP) //gcc-15.0.1 segfault here
{
};

template <class T> void Orbital_IBS_Common1<T>::Set(const LAParams& lap)
{
    itsLAParams=lap;
} 

template <class T>  LASolver<double>* Orbital_IBS_Common1<T>::CreateSolver() const
{
    LASolver<double>* las=LASolver<double>::Factory(itsLAParams);
    las->SetBasisOverlap(this->Overlap());
    return las;
}



template class TIBS_Common1<double>;
template class Orbital_IBS_Common1<double>;


// File: TBasisSetImplementation.C
module;
#include <tuple>
#include <iostream>
#include <vector>
import qchem.LAParams;
import qchem.LASolver;

module qchem.BasisSet.Internal.IrrepBasisSet;
import oml;

LAParams DefaultLAP({qchem::Lapack,qchem::SVD,1e-10,1e-12});
//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> Orbital_IBS_Common<T>::Orbital_IBS_Common()
    : itsLAParams      (DefaultLAP) //gcc-15.0.1 segfault here
{
};

template <class T> void Orbital_IBS_Common<T>::Set(const LAParams& lap)
{
    itsLAParams=lap;
} 

template <class T>  LASolver<double>* Orbital_IBS_Common<T>::CreateSolver() const
{
    LASolver<double>* las=LASolver<double>::Factory(itsLAParams);
    las->SetBasisOverlap(this->Overlap());
    return las;
}



template class IrrepBasisSet_Common<double>;
template class Orbital_IBS_Common<double>;


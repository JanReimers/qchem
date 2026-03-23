// File: TBasisSetImplementation.C
module;
#include <tuple>
#include <iostream>
#include <vector>
#include <memory>
#include <cassert>

import qchem.LAParams;
import qchem.LASolver_blaze;

module qchem.BasisSet.Internal.IrrepBasisSet;
import oml;

template <class T> IrrepBasisSet_Common<T>::sym_t IrrepBasisSet_Common<T>::GetSymmetry() const
    {
        assert(itsSymmetry);
        return itsSymmetry;
    }


LAParams DefaultLAP({qchem::Cholsky,1e-12});
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

template <class T>  LASolver_blaze<T>* Orbital_IBS_Common<T>::CreateSolver_blaze() const
{
    LASolver_blaze<T>* las=LASolver_blaze<T>::Factory(itsLAParams.BasisOrthoAlgorithm,itsLAParams.TruncationTolerance);
    las->SetBasisOverlap(LASolver_blaze<T>::convert(this->Overlap()));
    // std::cout << "Minimum singular value for basis set overlap= " << Min(las->Get_BS_Diagonal()) << std::endl;
    return las;
}



template class IrrepBasisSet_Common<double>;
template class Orbital_IBS_Common<double>;


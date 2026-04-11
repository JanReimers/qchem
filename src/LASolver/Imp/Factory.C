// File: LASolver/Imp/Factory.C
module;
#include <cassert>
module qchem.LASolver;
import qchem.LASolver.Internal.Lapack;

template <class T> LASolver<T>* LASolver<T>::
    Factory(qchem::Ortho ortho, double TruncationTolerance)
{
    LASolver<T>* ret=0;
    switch (ortho)
    {
    case qchem::Cholsky :
        ret=new LASolverCholsky<T>(TruncationTolerance);
        break;
    case qchem::Eigen :
        ret=new LASolverEigen<T>(TruncationTolerance);
        break;
    case qchem::SVD :
        ret=new LASolverSVD<T>(TruncationTolerance);
        break;
    break;
    }
    assert(ret);
    return ret;
}

template LASolver<double>* LASolver<double>::Factory(qchem::Ortho ortho, double TruncationTolerance);

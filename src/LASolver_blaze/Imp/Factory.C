// File: LASolver_blaze/Imp/Factory.C
module;
#include <cassert>
module qchem.LASolver;
import qchem.LASolver_blaze.Internal.Lapack;

template <class T> LASolver<T>* LASolver<T>::
    Factory(qchem::Ortho ortho, double TruncationTolerance)
{
    LASolver<T>* ret=0;
    switch (ortho)
    {
    case qchem::Cholsky :
        ret=new LASolverLapackCholsky_blaze<T>(TruncationTolerance);
        break;
    case qchem::Eigen :
        ret=new LASolverLapackEigen_blaze<T>(TruncationTolerance);
        break;
    case qchem::SVD :
        ret=new LASolverLapackSVD_blaze<T>(TruncationTolerance);
        break;
    break;
    }
    assert(ret);
    return ret;
}

template LASolver<double>* LASolver<double>::Factory(qchem::Ortho ortho, double TruncationTolerance);

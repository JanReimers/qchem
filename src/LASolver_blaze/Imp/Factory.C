// File: LASolver_blaze/Imp/Factory.C
module;
#include <cassert>
module qchem.LASolver_blaze;
import qchem.LASolver_blaze.Internal.Lapack;

template <class T> LASolver_blaze<T>* LASolver_blaze<T>::
    Factory(qchem::Ortho ortho, double TruncationTolerance)
{
    LASolver_blaze<T>* ret=0;
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

template LASolver_blaze<double>* LASolver_blaze<double>::Factory(qchem::Ortho ortho, double TruncationTolerance);

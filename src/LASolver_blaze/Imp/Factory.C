// File: Imp/Factory.C  
module;
#include <cassert>
module qchem.LASolver;
import qchem.LASolver.Internal.Lapack;
import qchem.LASolver.Internal.OML;

template <class T> LASolver<T>* LASolver<T>::
    Factory(const LAParams& lap)
{
    LASolver<T>* ret=0;
    switch (lap.LinearAlgebraPackage)
    {
    case qchem::OML:
        switch (lap.BasisOrthoAlgorithm)
        {
        case qchem::Cholsky :
            ret=new LASolverOMLCholsky<T>(lap);
            break;
        case qchem::Eigen :
            ret=new LASolverOMLEigen<T>(lap);
            break;
        case qchem::SVD :
            ret=new LASolverOMLSVD<T>(lap);
            break;
        }
        break;
   case qchem::Lapack:
        switch (lap.BasisOrthoAlgorithm)
        {
        case qchem::Cholsky :
            ret=new LASolverLapackCholsky<T>(lap);
            break;
        case qchem::Eigen :
            ret=new LASolverLapackEigen<T>(lap);
            break;
        case qchem::SVD :
            ret=new LASolverLapackSVD<T>(lap);
            break;
        break;
        }
    }
    assert(ret);
    return ret;
}

template LASolver<double>* LASolver<double>::Factory(const LAParams& lap);

// File: LASolver/Imp/Factory.C
module;
#include <cassert>
#include <iostream>
module qchem.LASolver;
import qchem.LASolver.Internal.Lapack;

using qchem::LASolverEigen;
using qchem::LASolverSVD;
using qchem::LASolverCholesky;
using qchem::LASolverCholeskyPivoted;
using qchem::LASolverAuto;

bool& qchem::ReportOverlapConditioning() { static bool on = false; return on; }

template <class T> LASolver<T>* LASolver<T>::
    Factory(qchem::Ortho ortho, double TruncationTolerance)
{
    LASolver<T>* ret = nullptr;
    switch (ortho)
    {
    case qchem::Cholesky:
        if (TruncationTolerance != 0)
        {
            std::cerr << "Warning: LASolverCholesky ignores TruncationTolerance ("
                      << TruncationTolerance << "); Cholesky does not truncate.\n";
        }
        ret = new LASolverCholesky<T>();
        break;
    case qchem::Eigen:   ret = new LASolverEigen<T>(TruncationTolerance); break;
    case qchem::SVD:     ret = new LASolverSVD  <T>(TruncationTolerance); break;
    case qchem::CholeskyPivoted: ret = new LASolverCholeskyPivoted<T>(TruncationTolerance); break;
    case qchem::Auto:
        if (TruncationTolerance != 0)
            std::cerr << "Warning: LASolverAuto ignores TruncationTolerance (" << TruncationTolerance
                      << "); Auto derives its own (gap detection).\n";
        ret = new LASolverAuto<T>();
        break;
    }
    assert(ret);
    return ret;
}

template LASolver<double>* LASolver<double>::Factory(qchem::Ortho, double);
template LASolver<dcmplx>* LASolver<dcmplx>::Factory(qchem::Ortho, double);

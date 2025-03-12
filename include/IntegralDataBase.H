// File: IntegralDataBase.H  Interface for a database of integrals over basis functions.
#ifndef _IntegralDataBase_H_
#define _IntegralDataBase_H_



#include "types.H"
#include "Imp/Containers/ERI4.H"

class IrrepIEClient;

//-----------------------------------------------------------------------------
//
//  Once supplied with an integral engine, this class manages storage
//  of all matrix elements.  All calculations are done by the integral
//  engine.  A smart implementation will request a calculation of the
//  data does not yet exist in the data base.
//
//
//  only build functions the have UT coverage.
//
// #define UT_COVERAGE_ONLY

// template <class T> class IntegralDataBase
// {
// public:
//     typedef Matrix<T>      Mat;
//     typedef SMatrix<T>     SMat;
//     typedef Vector<T>      Vec;
//     typedef Vector<double> RVec;
//     typedef std::vector<SMat> ERI3;
    
//     typedef const ScalarFunction<double> Rf;
//     typedef const ScalarFunction<T>      Sf;
//     typedef const VectorFunction<T>      Vf;
    
//     typedef const TIrrepBasisSet<T> bs_t;
//     typedef const IrrepIEClient iec_t;
    
//     virtual ~IntegralDataBase() {};

// };



#endif //_IntegralDataBase_H_

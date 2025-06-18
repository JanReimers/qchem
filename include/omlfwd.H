// File: oml/fwd.H  forward declare some oml classes.
#ifndef _oml_fwd_H_
#define _oml_fwd_H_

#include "oml/imp/index_t.h"

template <class T> class Matrix;
template <class T> class SMatrix;
template <class T> class DiagonalMatrix;
template <class T> class Vector;
template <class T> class Vector3D;

typedef Vector3D<double> RVec3; //Not really a fwd decalre!

#endif //_oml_fwd_H_
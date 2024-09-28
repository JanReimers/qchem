#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

#include "oml/imp/index_t.h"

//
//  Forward declarations for OML types.
//
template <class T> class Matrix;
template <class T> class SMatrix;
template <class T> class Vector;
template <class T> class Vector3D;

typedef Vector3D<double> RVec3;
//
//  Forward declarations for misc. types.
//
template <class T> class rc_ptr;

//
//  Forward declarations for functions types.
//
template <class T> class ScalarFunction;
template <class T> class VectorFunction;
//
//  Forward declarations for qchem abstract base classes.
//
class Atom;
class AngularMesh;
class BasisSet;
class BasisGroup;
class ChargeDensity;
class Cluster;
class ElectronDumper;
class EnergyLevel;
class ERIList;
class FittedCD;
class FittedFunction;
class FittedFunctionClient;
class Hamiltonian;
class HamiltonianTerm;
class Mesh;
class Orbital;
class OrbitalGroup;
class PlottableScalarFunction;
class PlottableVectorFunction;
class Plotter;
class QuantumNumber;
class RadialFunction;
class RadialMesh;
class SCFIterator;
class Spin;
class TotalEnergy;
class WaveFunction;


template <class T> class EigenSolver;
template <class T> class IntegralDataBase;
template <class T> class MatrixList;
template <class T> class TBasisFunction;
template <class T> class TBasisSet;
template <class T> class TOrbital;


#endif // TYPES_H_INCLUDED

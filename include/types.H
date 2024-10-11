#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

#include "oml/imp/index_t.h"

//
//  Forward declarations for OML types.
//
template <class T> class Matrix;
template <class T> class SMatrix;
template <class T> class DiagonalMatrix;
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
class IrrepBasisSet;
class BasisSet;
class ChargeDensity;
class Cluster;
class ElectronDumper;
class EnergyLevel;
class ERIList;
class ERI4;
class FittedCD;
class FittedFunction;
class FittedFunctionClient;
class Hamiltonian;
class HamiltonianTerm;
class IrrepIEClient;
class Mesh;
class Orbital;
class OrbitalGroup;
class PlottableScalarFunction;
class PlottableVectorFunction;
class Plotter;
class QuantumNumber;
//class RadialFunction;
class RadialMesh;
class SCFIterator;
class Spin;
class TotalEnergy;
class WaveFunction;

struct LinearAlgebraParams;
struct SCFIterationParams;


template <class T> class LASolver;
template <class T> class AnalyticIE;
template <class T> class IntegralDataBase;
template <class T> class MatrixList;
template <class T> class TBasisFunction;
template <class T> class TIrrepBasisSet;
template <class T> class TOrbital;


#endif // TYPES_H_INCLUDED

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
//  Forward declarations for functions types.
//
template <class T> class ScalarFunction;
template <class T> class VectorFunction;
//
//  Forward declarations for qchem abstract base classes.
//
class Atom;
class IrrepBasisSet;
class BasisSet;
class ChargeDensity;
class Cluster;
class ElectronConfiguration;
class EnergyLevel;
class EnergyLevels;
class ERI4;
class M4;
class ExFunctional;
class FittedCD;
class FittedFunction;
class FittedFunctionClient;
class Hamiltonian;
class HamiltonianTerm;
class Molecule;
class Mesh;
class MeshParams;
class Orbital;
class Orbitals;
class QuantumNumber;
class RadialMesh;
class SCFIterator;
class Spin;
class TotalEnergy;
class WaveFunction;

struct LAParams;
struct SCFIterationParams;


template <class T> class LASolver;
template <class T> class AnalyticIE;
template <class T> class IntegralDataBase;
template <class T> class TBasisFunction;
template <class T> class TIrrepBasisSet;
template <class T> class Orbital_IBS;
template <class T> class TOrbital;

namespace qchem
{
    // Integral types.
    enum IType3C {Overlap3C, Repulsion3C}; // <ab|c> and <ar|1/r12|c>
    enum IType2C {Overlap2C, Repulsion2C,Kinetic,Nuclear,RestMass, InvOverlap, InvRepulsion, Charge, Normalization};
    enum IType   {Overlap1, Kinetic1, Nuclear1, RestMass1,Repulsion1,InvOverlap1, InvRepulsion1, Charge1};
}

#endif // TYPES_H_INCLUDED

#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

#include "oml/imp/index_t.h"

/*! \mainpage qchem .
 *
 * \section overview_sec Overview 
 * 
 * This code is an experiment in using SOLID principles of object oriented design
 * for a quantum chemistry application.  For example a useful and obvious abstraction is a
 * basis set, which can be implmented as polarized gaussians, spherical gaussiains, Slater functions,
 * B-plines, etc. 
 *
 * \section ac_sec Acronyms
 * 
 * In order to avoid very long class names we will use a number of acronyms and short forms. These acronyms can
 * also be used in file names.
 * 
 * Acronym  | Expansion         |  Comments |
 * ---------| ----------------- | --------------------------- |
 * BF       | Basis Function    | Gaussian, Slater, etc. |
 * BS       | Basis Set         | Full basis set == list of IBSs. |
 * IBS      | Irrep Basis Set   | List of Basis Functions. |
 * FBS      | Fit Basis Sets    | For DFT where we fit Vxc and CD |
 * IE       | Integral Engine   | Calculates tables of integrals. |
 * DB       | DataBase          | Caches integral tables. |
 * CD       | Charge Density    | |
 * DM       | Density Matrix    | |
 * HF       | Hartree-Fock      | |
 * DHF      | Dirac-Hartree-Fock      | |
 * DFT      | Density Functional Theory | |
 * WF       | Wave Function     | |
 * HT       | Hamiltonian Term  | i.e. Kinetic, Ven, Vee, Vxc etc |
 * SG       | Spherical Gaussian | spdfg gaussian functions used for atom calculations.
 * SS       | Spherical Slater   | spdfg Slater functions used for atom calculations.
 * PG       | Polarized Gaussian | Polarized (x,y,z,xx,xy,xz,..) Gaussians used for molecular calculations
 * QN       | Quantum Number    | Usually a group of QNs like {n,l,ml,ms}
 * EL       | Energy level      | Eigen energies for orbitals.
 * EC       | Electron Configuration | Specifies occupied orbitals, i.e. [Xe]6s^2 4f^1 5d^1 |
 * ERI      | Electron Repulsion Integral   | |
 * Occ      | Occupation        | |
 * RKB      | Restricted Kinetic Balance | Used for relativistic Basis Sets |
 * 
 * 
 * Irrep =  Irreducible Representation
*/


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
//  Some stuff in std::
//
namespace std 
{
    template <typename T> class valarray;
}
//
//  Forward declarations for function types.
//
template <class T> class ScalarFunction;
template <class T> class VectorFunction;
//
//  Forward declarations for qchem abstract base classes.
//
class Atom;
class IrrepBasisSet;
class BasisFunction;
class BasisSet;
class Cluster;
class Dynamic_HT;
class ElectronConfiguration;
class EnergyLevel;
class EnergyLevels;
class ERI4;
class DM_CD;
class M4;
class ExFunctional;
class Fit_IBS;
class FittedCD;
class FittedFunction;
class FittedFunctionClient;
class Hamiltonian;
class Static_HT;
class Irrep_QNs;
class Irrep_WF;
class Molecule;
class Mesh;
class MeshParams;
class Orbital_IBS;
class Orbital;
class Orbitals;
class Orbital_QNs;
class SCFIrrepAccelerator;
class SCFAccelerator;
class Symmetry;
enum class Spin;
class RadialMesh;
class SCFIterator;
class EnergyBreakdown;
class WaveFunction;

struct LAParams;
struct SCFIterationParams;


template <class T> class LASolver;
template <class T> class AnalyticIE;
template <class T> class TBasisFunction;
template <class T> class TIrrepBasisSet;
template <class T> class TOrbital_IBS;
template <class T> class TOrbital_DFT_IBS;
template <class T> class TOrbital_HF_IBS;
template <class T> class Orbital_RKBL_IBS;
template <class T> class Orbital_RKBS_IBS;
template <class T> class TOrbital;

namespace qchem
{
    // Integral types.
    enum IType3C {Overlap3C, Repulsion3C}; // <ab|c> and <ar|1/r12|c>
    enum IType2C {Overlap2C, Repulsion2C,Grad2,Nuclear,RestMass, InvOverlap, InvRepulsion, Charge, Normalization,NumCharge,NumNormalization,NumOverlap};
    enum IType   {Overlap1, Grad2_1, Nuclear1, RestMass1,Repulsion1,InvOverlap1, InvRepulsion1, Charge1};
}

#endif // TYPES_H_INCLUDED

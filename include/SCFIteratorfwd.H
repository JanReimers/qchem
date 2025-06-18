#ifndef TYPES_H_INCLUDED
#define TYPES_H_INCLUDED

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

class SCFIterator;
struct SCFParams;

#endif // TYPES_H_INCLUDED

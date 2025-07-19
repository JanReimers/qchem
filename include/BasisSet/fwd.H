// File: BasisSet/fwd.H  forward declare some BasisSet classes.
#ifndef _BasisSet_fwd_H_
#define _BasisSet_fwd_H_

//Super matrix containers
class ERI4;
class M4;

class IrrepBasisSet;
class Fit_IBS;
class Orbital_IBS;

class BasisFunction;

template <class T> class TIrrepBasisSet;
template <class T> class     TOrbital_IBS;
template <class T> class TOrbital_DFT_IBS;
template <class T> class  TOrbital_HF_IBS;
template <class T> class Orbital_RKBL_IBS;
template <class T> class Orbital_RKBS_IBS;

template <class T> class TBasisFunction;

namespace qchem
{
    // Integral types.
    enum IType3C {Overlap3C, Repulsion3C}; // <ab|c> and <ar|1/r12|c>
    enum IType2C {Overlap2C, Repulsion2C,Grad2,Nuclear,RestMass, InvOverlap, InvRepulsion, Charge, Normalization,NumCharge,NumNormalization,NumOverlap};
    enum IType   {Overlap1, Grad2_1, Nuclear1, RestMass1,Repulsion1,InvOverlap1, InvRepulsion1, Charge1};
}

#endif //_BasisSet_fwd_H_
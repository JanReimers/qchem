// File: BasisSet/Imp/IntegralEnums.C  Define some enums for integral types;
export module qchem.BasisSet.Internal.IntegralEnums;

export namespace qchem
{
    // Integral types.
    enum IType3C {Overlap3C, Repulsion3C}; // <ab|c> and <ar|1/r12|c>
    enum IType2C {Overlap2C, Repulsion2C,Grad2,Nuclear,RestMass, InvOverlap, InvRepulsion, Charge, Normalization,NumCharge,NumNormalization,NumOverlap};
    enum IType   {Overlap1, Grad2_1, Nuclear1, RestMass1,Repulsion1,InvOverlap1, InvRepulsion1, Charge1};
}


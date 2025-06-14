// File: SphericalSymmetryQN.N  A Quantum Number for atomic (spherical) symmetry



#include "Symmetry/Yl.H"
#include "Imp/BasisSet/Atom/EC.H"
#include <iostream>
#include <cassert>

Yl_Sym::Yl_Sym()
    : itsL(0)
{};

Yl_Sym::Yl_Sym(int theL)
    : itsL(theL)
{};

 size_t Yl_Sym::SequenceIndex() const //Used for op<
 {
    return itsL;
 }


int Yl_Sym::GetDegeneracy() const
{
    return 2*itsL+1;
}

ElCounts_l Yl_Sym::GetN(const ElCounts& ec) const
{
    return ElCounts_l{ec.N[itsL],ec.Nu[itsL]}; // Should be total core+valance 
}


std::string SPDFG[]={"s","p","d","f","g"};

std::ostream& Yl_Sym::Write(std::ostream& os) const
{
    return os << SPDFG[itsL] << " ";
}


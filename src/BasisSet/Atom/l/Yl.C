// File: SphericalSymmetryQN.N  A Quantum Number for atomic (spherical) symmetry



#include "Imp/BasisSet/Atom/l/Yl.H"
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

bool Yl_Sym::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const Yl_Sym*>(&b)!=0;
}

bool Yl_Sym::Match(const Symmetry& qn) const
{
    const Yl_Sym* aqn = dynamic_cast<const Yl_Sym*>(&qn);
    assert(aqn);
    return itsL==aqn->itsL;
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

Angular_Sym* Yl_Sym::Clone() const
{
    return new Yl_Sym(*this);
}


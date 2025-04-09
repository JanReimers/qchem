// File: SphericalSymmetryQN.N  A Quantum Number for atomic (spherical) symmetry



#include "Imp/Symmetry/YlQN.H"
#include "Imp/WaveFunction/Atom_EC.H"
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
    int nl=ec.N[itsL];
    if (ec.Nv[itsL]==0) return ElCounts_l{nl,0};//std::make_pair(nl,0);
    assert(nl!=0);
    // Handle partial shells
    int nlu=1; //# unpaired in shell l. 
    if (itsL==1) // p is partial.
    {
        assert(ec.Nv[2]==0); //No partial D orbital
        nlu=ec.NUnpaired-ec.Nv[0];
    }
    else if (itsL==2) // d is partial.
    {
        assert(ec.Nv[1]==0); //p better be full
        if (ec.Nv[itsL]>1) nlu=ec.NUnpaired-ec.Nv[0];            
    }
    else if(itsL==3) // f is partial.
    {
        
        assert(ec.Nv[0]==0); //If f is Partial s must be full.
        assert(ec.Nv[1]==0); //If f is Partial p must be full.
        nlu=ec.NUnpaired-ec.Nv[2];
        assert(nlu>=0);
    }
    return ElCounts_l{nl,nlu};//std::make_pair(nl,nlu);
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


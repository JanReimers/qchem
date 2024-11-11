// File: SphericalSymmetryQN.N  A Quantum Number for atomic (spherical) symmetry



#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>

SphericalSymmetryQN::SphericalSymmetryQN()
    : itsL(0)
{};

SphericalSymmetryQN::SphericalSymmetryQN(int theL)
    : itsL(theL)
{};

bool SphericalSymmetryQN::Match(const QuantumNumber& qn) const
{
    const SphericalSymmetryQN* aqn = dynamic_cast<const SphericalSymmetryQN*>(&qn);
    assert(aqn);
    return itsL==aqn->itsL;
}

int SphericalSymmetryQN::GetDegeneracy() const
{
    return 2*itsL+1;
}

std::pair<int,int> SphericalSymmetryQN::GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const
{
    int nl=N[itsL];
    if (Nv[itsL]==0) return std::make_pair(nl,0);
    assert(nl!=0);
    // Handle partial shells
    int nlu=1; //# unpaired in shell l. 
    if (itsL==1) // p is partial.
    {
        assert(Nv[2]==0); //No partial D orbital
        nlu=NUnpaired-Nv[0];
    }
    else if (itsL==2) // d is partial.
    {
        assert(Nv[1]==0); //p better be full
        if (Nv[itsL]>1) nlu=NUnpaired-Nv[0];            
    }
    else if(itsL==3) // f is partial.
    {
        
        assert(Nv[0]==0); //If f is Partial s must be full.
        assert(Nv[1]==0); //If f is Partial p must be full.
        nlu=NUnpaired-Nv[2];
        assert(nlu>=0);
    }
    return std::make_pair(nl,nlu);
}


std::string SPDFG[]={"s","p","d","f","g"};

std::ostream& SphericalSymmetryQN::Write(std::ostream& os) const
{
    UniqueID::Write(os);
    if (StreamableObject::Binary())
        BinaryWrite(itsL,os);
    if (StreamableObject::Ascii())
        os << itsL << " ";
    if (StreamableObject::Pretty())
        os << SPDFG[itsL] << " ";
    return os;
}

std::istream& SphericalSymmetryQN::Read (std::istream& is)
{
    UniqueID::Read(is);
    if (StreamableObject::Binary())
        BinaryRead(itsL,is);
    else
    {
        is >> itsL;
    }
    return is;
}

QuantumNumber* SphericalSymmetryQN::Clone() const
{
    return new SphericalSymmetryQN(*this);
}


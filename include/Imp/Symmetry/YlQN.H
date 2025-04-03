// File: Symmetry/YlmQN.H  A Quantum Number for atomic (spherical) symmetry
#ifndef _Symmetry_YlQN_H_
#define _Symmetry_YlQN_H_

#include "Imp/Symmetry/AngularQN.H"

//---------------------------------------------------------------------------------
//
//  The relevant quantum number for spherical symmetry is the magnitude of the
//  orbital angular momentum L.
//

class YlQN
    : public virtual AngularQN
{
public:
    YlQN(        );
    YlQN(int theL);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual bool MatchType(const Symmetry&) const;
    virtual bool Match(const Symmetry&) const;
    virtual int GetDegeneracy() const;
    virtual Symmetry* AddPrincipleQN(int index) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual AngularQN* Clone(        ) const;

    virtual int     GetL() const {return itsL;}
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

protected:
    int itsL;
    static const int LMax=3;
};

#endif //_Symmetry_YlQN_H_

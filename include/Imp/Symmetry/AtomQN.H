// File: Symmetry/AtomQN.H  Spherical harmonic Ylm symmetry plus princible (radial QN).
#ifndef _Symmetry_AtomQN_H_
#define _Symmetry_AtomQN_H_

#include <Symmetry.H>
#include <tuple>

class AngularQN;

//---------------------------------------------------------------------------------
//
//  Atomic QN is a combination of the principle QN and the angular QN.
//

class AtomQN
    : public virtual Symmetry
{
public:
    AtomQN(            );
    AtomQN(int n, const AngularQN&);
    AtomQN(const AtomQN&);
    virtual ~AtomQN();
    
    virtual size_t SequenceIndex() const; //Used for op<
    virtual bool MatchType(const Symmetry&) const;
    
    virtual bool Match(const Symmetry&) const;
    virtual int GetDegeneracy() const;
    virtual Symmetry* AddPrincipleQN(int index) const;
 
    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual Symmetry* Clone(        ) const;
private:
    int l,n;
    AngularQN* itsAngularQN; //Yl Ylml, or relativistic Ok, Okmj
};

#endif //_Symmetry_AtomQN_H_

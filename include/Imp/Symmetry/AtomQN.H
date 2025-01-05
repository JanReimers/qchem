// File: Symmetry/AtomQN.H  Spherical harmonic Ylm symmetry plus princible (radial QN).
#ifndef _Symmetry_AtomQN_H_
#define _Symmetry_AtomQN_H_

#include <QuantumNumber.H>
#include <tuple>

class AngularQN;

//---------------------------------------------------------------------------------
//
//  Atomic QN is a combination of the principle QN and the angular QN.
//

class AtomQN
    : public virtual QuantumNumber
{
public:
    AtomQN(            );
    AtomQN(int n, const AngularQN&);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;
    virtual QuantumNumber* AddPrincipleQN(int index) const;
 
    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;
private:
    int l,n;
    AngularQN* itsAngularQN; //Yl Ylml, or relativistic Ok, Okmj
};

#endif //_Symmetry_AtomQN_H_

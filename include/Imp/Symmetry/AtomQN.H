// File: Symmetry/AtomQN.H  Spherical harmonic Ylm symmetry plus princible (radial QN).
#ifndef _Symmetry_AtomQN_H_
#define _Symmetry_AtomQN_H_

#include <QuantumNumber.H>

class YlQN;
class Omega_kmjQN;
//---------------------------------------------------------------------------------
//
//  angular momentum L and m.
//

class AtomQN
    : public virtual QuantumNumber
{
public:
    AtomQN(            );
    AtomQN(int n, const YlQN&);
    AtomQN(int n, const Omega_kmjQN&);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;
    virtual QuantumNumber* AddPrincipleQN(int index) const;
    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;
private:
    int l,n;
    QuantumNumber* itsAngularQN; //Yl Ylml, or relativistic Ok, Okmj
};

#endif //_Symmetry_AtomQN_H_

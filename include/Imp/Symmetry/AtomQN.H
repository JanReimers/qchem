// File: Symmetry/AtomQN.H  Spherical harmonic Ylm symmetry plus princible (radial QN).
#ifndef _Symmetry_AtomQN_H_
#define _Symmetry_AtomQN_H_

#include "Imp/Symmetry/YlmQN.H"


//---------------------------------------------------------------------------------
//
//  angular momentum L and m.
//

class AtomQN
    : public virtual QuantumNumber
    , public YlmQN
{
public:
    AtomQN(            );
    AtomQN(int n, const YlmQN&);
    AtomQN(int n, int l, int m);

    virtual bool Match(const QuantumNumber&) const;
    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;
private:
    int n;
};

#endif //_Symmetry_AtomQN_H_

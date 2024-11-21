// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/Symmetry/AtomQN.H"
#include <iostream>
#include <cassert>

using std::cout;
using std::endl;

AtomQN::AtomQN(): YlmQN(0,0),n(0) {};

AtomQN::AtomQN(int _n, int _l, int _m) : YlmQN(_l,_m), n(_n) {};

AtomQN::AtomQN(int _n, const YlmQN& Ylm) : YlmQN(Ylm), n(_n+itsL) {};
AtomQN::AtomQN(int _n, const YlQN& Yl) : YlmQN(Yl.GetL(),0), n(_n+itsL) {};


bool AtomQN::Match(const QuantumNumber& qn) const
{
    const AtomQN* yqn = dynamic_cast<const AtomQN*>(&qn);
    assert(yqn);
    return n==yqn->n && itsL==yqn->itsL && m==yqn->m;
}



std::ostream& AtomQN::Write(std::ostream& os) const
{
    os << n << " ";
    return YlmQN::Write(os);
}

std::istream& AtomQN::Read (std::istream& is)
{
    return is;
}

QuantumNumber* AtomQN::Clone() const
{
    return new AtomQN(*this);
}


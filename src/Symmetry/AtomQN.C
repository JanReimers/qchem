// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/Symmetry/AtomQN.H"
#include "Imp/Symmetry/YlmQN.H"
#include "Imp/Symmetry/OkmjQN.H"
#include <iostream>
#include <cassert>

using std::cout;
using std::endl;

AtomQN::AtomQN(): l(0),n(0),itsAngularQN(0) {};

// AtomQN::AtomQN(int _n, const YlQN& Yl) 
// : l(Yl.GetL()), n(_n+l), itsAngularQN(Yl.Clone()) {};

// AtomQN::AtomQN(int _n, const Omega_kmjQN& Okmj)
// : l(Okmj.GetL()), n(_n+l), itsAngularQN(Okmj.Clone()) {};

AtomQN::AtomQN(int _n, const AngularQN& aqn)
: l(aqn.GetL()), n(_n+l), itsAngularQN(aqn.Clone()) {};

bool AtomQN::Match(const QNs& qn) const
{
    const AtomQN* yqn = dynamic_cast<const AtomQN*>(&qn);
    assert(yqn);
    return n==yqn->n && l == yqn->l;
}

int AtomQN::GetDegeneracy() const
{
    return itsAngularQN->GetDegeneracy();
}

QNs* AtomQN::AddPrincipleQN(int index) const
{
    return itsAngularQN->AddPrincipleQN(index);
}


std::ostream& AtomQN::Write(std::ostream& os) const
{
    return os << n << *itsAngularQN;
}

std::istream& AtomQN::Read (std::istream& is)
{
    return is;
}

QNs* AtomQN::Clone() const
{
    return new AtomQN(*this);
}


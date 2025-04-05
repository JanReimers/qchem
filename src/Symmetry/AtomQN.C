// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.



#include "Imp/Symmetry/AtomQN.H"
#include "Imp/Symmetry/YlmQN.H"
#include "Imp/Symmetry/OkmjQN.H"
#include <iostream>
#include <cassert>

using std::cout;
using std::endl;
const int n_max=20;

AtomQN::AtomQN(): l(0),n(0),itsAngularQN(0) {};

// AtomQN::AtomQN(int _n, const YlQN& Yl) 
// : l(Yl.GetL()), n(_n+l), itsAngularQN(Yl.Clone()) {};

// AtomQN::AtomQN(int _n, const Omega_kmjQN& Okmj)
// : l(Okmj.GetL()), n(_n+l), itsAngularQN(Okmj.Clone()) {};

AtomQN::AtomQN(int _n, const AngularQN& aqn)
: l(aqn.GetL()), n(_n+l), itsAngularQN(aqn.Clone()) {};

AtomQN::AtomQN(const AtomQN& aqn)
: l(aqn.l), n(aqn.n), itsAngularQN(aqn.itsAngularQN->Clone()) {};

AtomQN::~AtomQN()
{
    if (itsAngularQN) delete itsAngularQN;
}

size_t AtomQN::SequenceIndex() const //Used for op<
 {
    assert(n<=n_max);
    return itsAngularQN->SequenceIndex()*n_max+(n-1);
 }

bool AtomQN::MatchType(const Symmetry& b) const
{
    return dynamic_cast<const AtomQN*>(&b)!=0;
}

bool AtomQN::Match(const Symmetry& qn) const
{
    const AtomQN* yqn = dynamic_cast<const AtomQN*>(&qn);
    assert(yqn);
    return n==yqn->n && l == yqn->l;
}

int AtomQN::GetDegeneracy() const
{
    return itsAngularQN->GetDegeneracy();
}

Symmetry* AtomQN::AddPrincipleQN(int index) const
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

Symmetry* AtomQN::Clone() const
{
    return new AtomQN(*this);
}


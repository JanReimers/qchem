// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.

#include <Orbital_QNs.H>
#include <QuantumNumber.H>
#include <cassert>
#include <iostream>

const size_t n_max=20;
const size_t ms_max=3; //three states Up/Down and None.


Orbital_QNs::Orbital_QNs(size_t _n, Spin _ms,const QNs* _sym)
: n(_n), ms(_ms), sym(_sym->Clone())
{
    assert(sym);
    assert(n<=n_max);
}

Orbital_QNs::~Orbital_QNs()
{
    delete sym;
}

size_t Orbital_QNs::SequenceIndex() const
{
    assert(n-1<n_max);
    assert(ms.SequenceIndex()<ms_max);
    size_t is=1;//sym->SequenceIndex();
    return (is*n_max+(n-1))*ms_max+ms.SequenceIndex();
}

bool Orbital_QNs::Match(const Orbital_QNs& b) const
{
    return n==b.n && ms==b.ms && (*sym==*b.sym);
}
bool Orbital_QNs::MatchType(const Orbital_QNs&) const
{
    return false;//sym->MatchType(*b.sym);
}
    
bool operator<(const Orbital_QNs& a, const Orbital_QNs& b)
{
    return a.SequenceIndex()<b.SequenceIndex();
}
    
size_t Orbital_QNs::GetDegeneracy() const
{
    return sym->GetDegeneracy()*ms.GetDegeneracy();
}

std::ostream& Orbital_QNs::Write(std::ostream& os) const
{
    return os << "n=" << n << " ms=" << ms << " sym=" << sym;
}
    
Orbital_QNs* Orbital_QNs::Clone() const
{
    return new Orbital_QNs(n,ms,sym);
}
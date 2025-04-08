// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.

#include <Orbital_QNs.H>
#include <Symmetry.H>
#include <cassert>
#include <iostream>

const size_t n_max=300; //Max principle QN.
const size_t ms_max=3; //three states Up/Down and None.

Irrep_QNs::Irrep_QNs(Spin _ms,const Symmetry* _sym) 
    : ms(_ms)
    , sym(_sym->Clone()) 
{
    assert(sym);
}

Irrep_QNs::Irrep_QNs(const Irrep_QNs& qns)
    : ms(qns.ms)
    , sym(qns.sym ? qns.sym->Clone() : 0)
{
    assert(sym);
}

Irrep_QNs::~Irrep_QNs()
{
    assert(sym);
    if (sym) delete sym;
}
bool Irrep_QNs::Match(const Irrep_QNs& b) const
{
    return ms==b.ms && (*sym==*b.sym);
}
size_t Irrep_QNs::SequenceIndex() const
{
    assert(ms.SequenceIndex()<ms_max);
    size_t is=sym->SequenceIndex();
    return is*ms_max+ms.SequenceIndex();
}
size_t Irrep_QNs::GetDegeneracy() const
{
    return sym->GetDegeneracy()*ms.GetDegeneracy();
}
std::ostream& Irrep_QNs::Write(std::ostream& os) const
{
    return os << " ms=" << ms << " sym=" << *sym;
}

Orbital_QNs::Orbital_QNs(size_t _n, Spin _ms,const Symmetry* _sym)
: n(_n), ms(_ms), sym(_sym->Clone())
{
    assert(sym);
    assert(n<=n_max);
}
Orbital_QNs::Orbital_QNs(const Orbital_QNs& qns)
: n(qns.n)
, ms(qns.ms)
, sym(qns.sym ? qns.sym->Clone() : 0)
{
    assert(sym);
}
Orbital_QNs::~Orbital_QNs()
{
    if (sym) delete sym;
}
size_t Orbital_QNs::SequenceIndex() const
{
    assert(n-1<n_max);
    assert(ms.SequenceIndex()<ms_max);
    Irrep_QNs iqns(ms,sym);
    return iqns.SequenceIndex()*n_max + (n-1);
}
bool Orbital_QNs::Match(const Orbital_QNs& b) const
{
    return n==b.n && ms==b.ms && (*sym==*b.sym);
}
bool Orbital_QNs::MatchNoSpin(const Orbital_QNs& b) const
{
    return n==b.n && (*sym==*b.sym);
}
size_t Orbital_QNs::GetDegeneracy() const
{
    return sym->GetDegeneracy()*ms.GetDegeneracy();
}

std::ostream& Orbital_QNs::Write(std::ostream& os) const
{
    return os << n << *sym;
}
    
// Orbital_QNs* Orbital_QNs::Clone() const
// {
//     return new Orbital_QNs(n,ms,sym);
// }
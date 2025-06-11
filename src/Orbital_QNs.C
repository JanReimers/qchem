// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.

#include <Orbital_QNs.H>
#include <Symmetry.H>
#include <cassert>
#include <iostream>


Orbital_QNs::Orbital_QNs(size_t _n, Spin _ms,const sym_t& _sym)
: n(_n), ms(_ms), sym(_sym)
{
    assert(sym);
    assert(n<=Irrep_QNs::n_max);
}

Orbital_QNs::Orbital_QNs(size_t _n, const Irrep_QNs& irr)
: n(_n)
, ms(irr.ms)
, sym(irr.sym)
{
    assert(sym);
}

Orbital_QNs::~Orbital_QNs()
{
}
size_t Orbital_QNs::SequenceIndex() const
{
    assert(n-1<Irrep_QNs::n_max);
    assert(::SequenceIndex(ms)<Irrep_QNs::ms_max);
    Irrep_QNs iqns(ms,sym);
    return iqns.SequenceIndex()*Irrep_QNs::n_max + (n-1);
}
size_t Orbital_QNs::GetDegeneracy() const
{
    return sym->GetDegeneracy()*::GetDegeneracy(ms);
}

std::ostream& Orbital_QNs::Write(std::ostream& os) const
{
    return os << n+sym->GetPrincipleOffset() << *sym;
}
    

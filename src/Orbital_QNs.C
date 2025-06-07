// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.

#include <Orbital_QNs.H>
#include <Symmetry.H>
#include <cassert>
#include <iostream>

const size_t n_max=300; //Max principle QN.
const size_t ms_max=3; //three states Up/Down and None.

Irrep_QNs::Irrep_QNs(Spin _ms,const sym_t& _sym) 
    : ms(_ms)
    , sym(_sym) 
{
    assert(sym);
}



Irrep_QNs::~Irrep_QNs()
{
   
}
size_t Irrep_QNs::SequenceIndex() const
{
    assert(::SequenceIndex(ms)<ms_max);
    size_t is=sym->SequenceIndex();
    return is*ms_max+::SequenceIndex(ms);
}
size_t Irrep_QNs::GetDegeneracy() const
{
    return sym->GetDegeneracy()*::GetDegeneracy(ms);
}
std::ostream& Irrep_QNs::Write(std::ostream& os) const
{
    return os << " ms=" << static_cast<int>(ms) << " sym=" << *sym;
}

Orbital_QNs::Orbital_QNs(size_t _n, Spin _ms,const sym_t& _sym)
: n(_n), ms(_ms), sym(_sym)
{
    assert(sym);
    assert(n<=n_max);
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
    assert(n-1<n_max);
    assert(::SequenceIndex(ms)<ms_max);
    Irrep_QNs iqns(ms,sym);
    return iqns.SequenceIndex()*n_max + (n-1);
}
size_t Orbital_QNs::GetDegeneracy() const
{
    return sym->GetDegeneracy()*::GetDegeneracy(ms);
}

std::ostream& Orbital_QNs::Write(std::ostream& os) const
{
    return os << n+sym->GetPrincipleOffset() << *sym;
}
    

module;
#include <cassert>
#include <iostream>

module qchem.Symmetry.Orbital;

const size_t Orbital_QNs::n_max=300; //Max principle QN.

Orbital_QNs::Orbital_QNs(size_t _n, Spin _ms,const sym_t& _sym)
: Irrep_QNs(_ms,_sym), n(_n)
{
    assert(n<=n_max);
}

Orbital_QNs::Orbital_QNs(size_t _n, const Irrep_QNs& irr)
: Irrep_QNs(irr), n(_n)
{}

Orbital_QNs::~Orbital_QNs()
{
}
size_t Orbital_QNs::SequenceIndex() const
{
    assert(n-1<n_max);
    return Irrep_QNs::SequenceIndex()*n_max + (n-1);
}

std::ostream& Orbital_QNs::Write(std::ostream& os) const
{
    os << n;
    return  Irrep_QNs::Write(os);
}
    

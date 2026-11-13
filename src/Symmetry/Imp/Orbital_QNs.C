module;
#include <cassert>
#include <iostream>

module qchem.Symmetry.Orbital;

namespace qchem {

Orbital_QNs::Orbital_QNs(size_t _n, Spin _ms,const sym_t& _sym)
: Irrep(_ms,_sym), n(_n)
{}

Orbital_QNs::Orbital_QNs(size_t _n, const Irrep& irr)
: Irrep(irr), n(_n)
{}

Orbital_QNs::~Orbital_QNs()
{
}
std::ostream& Orbital_QNs::Write(std::ostream& os) const
{
    os << n;
    return  Irrep::Write(os);
}
    

} // namespace qchem
// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.

#include <cassert>
#include <iostream>
#include <Common/pmstream.h>
#include <Symmetry/Irrep_QNs.H>
import qchem.Symmetry;

const size_t Irrep_QNs::ms_max=3; //three states Up/Down and None.

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

std::string spins[]={"↓"," ","↑"};
std::ostream& Irrep_QNs::Write(std::ostream& os) const
{
    return os  << (*sym) << spins[static_cast<int>(ms)];
}

    

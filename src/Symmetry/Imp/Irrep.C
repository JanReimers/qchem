// File: Irrep.C  Combine Symmetry with Spin.
module;
#include <cassert>
#include <iostream>
module qchem.Symmetry.Irrep;
import qchem.Streamable;
import qchem.Strings;

namespace qchem {

const size_t Irrep::ms_max=3; //three states Up/Down and None.

Irrep::Irrep(Spin _ms,const sym_t& _sym) 
    : ms(_ms)
    , sym(_sym) 
{
    assert(sym);
}



Irrep::~Irrep()
{
   
}
size_t Irrep::SequenceIndex() const
{
    assert(qchem::SequenceIndex(ms)<ms_max);
    size_t is=sym->SequenceIndex();
    return is*ms_max+qchem::SequenceIndex(ms);
}
size_t Irrep::GetDegeneracy() const
{
    return sym->GetDegeneracy()*qchem::GetDegeneracy(ms);
}

std::ostream& Irrep::Write(std::ostream& os) const
{
    return os  << (*sym) << spins[static_cast<int>(ms)];
}

    

} // namespace qchem
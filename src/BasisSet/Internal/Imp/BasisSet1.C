// File: BasisSet/Imp/BS_Common.C  Common implementation for BasisSet
module;
#include <cassert>
#include <memory>
import qchem.LAParams;

module qchem.BasisSet.Internal.Common1;
import qchem.stl_io;
import qchem.IrrepBasisSet;
import qchem.Streamable;

void BS_Common1::Insert(bs_t* bs)
{
    assert(bs);
    itsBasisSets.push_back(std::unique_ptr<bs_t>(bs));
}

size_t BS_Common1::GetNumFunctions() const
{
    size_t ret=0;
    for (auto& bs:*this) 
        ret+=bs->GetNumFunctions();
    return ret;
}

BasisSet1::irrepv_t BS_Common1::GetIrreps(const Spin& ms) const
{
    irrepv_t irrepv;
    for (auto& b:itsBasisSets) irrepv.push_back(b->GetIrrep(ms));
    return irrepv;
}
//
//  StreamableObject stuff.
//
std::ostream&  BS_Common1::Write(std::ostream&  os  ) const
{
    return os << itsBasisSets;
}


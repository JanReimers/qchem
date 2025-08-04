// File: BasisSet/Imp/BS_Common.C  Common implementation for BasisSet
module;
#include <cassert>
#include <memory>
import qchem.LAParams;

module qchem.BasisSet.Internal.Common;
import qchem.stl_io;
import qchem.IrrepBasisSet;

void BS_Common::Insert(bs_t* bs)
{
    assert(bs);
    itsBasisSets.push_back(std::unique_ptr<bs_t>(bs));
}

void BS_Common::Set(const LAParams& lap)
{
    for (auto& b:itsBasisSets) b->Set(lap);
}

size_t BS_Common::GetNumFunctions() const
{
    size_t ret=0;
    for (auto& bs:*this) 
        ret+=bs->GetNumFunctions();
    return ret;
}

BasisSet::symv_t BS_Common::GetSymmetries  () const
{
    symv_t symv;
    for (auto& b:itsBasisSets) symv.push_back(b->GetSymmetry());
    return symv;
}

//
//  StreamableObject stuff.
//
std::ostream&  BS_Common::Write(std::ostream&  os  ) const
{
    return os << itsBasisSets;
}


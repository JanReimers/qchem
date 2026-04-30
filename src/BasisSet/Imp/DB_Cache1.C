// File: BasisSet/Imp/DB_Cache1.C Global integrals cache allow data sharing between separate runs.
module;
#include <map>
#include <cassert>
#include <vector>
#include <string>
module qchem.BasisSet.DB_Cache1;
// import qchem.BasisSet.Internal.ERI4;
// import qchem.BasisSet.Internal.ERI3;
// import qchem.BasisSet.Internal.IntegralEnums;

size_t Cache41::Register(const std::string& bf_id)
{
    size_t index;
    auto i=itsUniqueBFs.find(bf_id);
    if (i==itsUniqueBFs.end())
        itsUniqueBFs[bf_id]=index=itsUniqueBFs.size();
    else
        index=i->second;
    return index;
}
// File: IEClient.H  Array like view of an Irrep basis set.
module;
#include <cstddef>
export module qchem.BasisSet.Imp.IEClient;
import Common.UniqueID; 

//
//  Integral DB and engines only sees this
//
export class IrrepIEClient //Client for and Irrep basis set.
: public virtual UniqueID
{
public:
    virtual ~IrrepIEClient() {};
    virtual size_t size() const=0;
};


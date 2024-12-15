
#include "Imp/BasisSet/Slater_mj/IEClient.H"

namespace Slater_mj
{
    
void Dirac_IrrepIEClient::Init(const Large_IrrepIEClient* liec,const Small_IrrepIEClient* siec)
{
    itsLargeIEC=liec;
    itsSmallIEC=siec;
    assert(itsLargeIEC);
    assert(itsSmallIEC);
}

size_t Dirac_IrrepIEClient::size() const
{
    assert(itsLargeIEC);
    assert(itsSmallIEC);
    return itsLargeIEC->size()+itsSmallIEC->size();
}

} //namespace

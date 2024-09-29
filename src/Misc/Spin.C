// File: Spin.C  A 3 state variable indicating the spin polarization.

#include "Misc/Spin.H"
#include "oml/imp/stream.h"
#include "oml/imp/binio.h"
#include <cassert>

std::ostream& Spin::Write(std::ostream& os) const
{
    if (StreamableObject::Binary())
    {
        int temp=static_cast<int>(itsState);
        BinaryWrite(temp,os);
    }
    else
        os << itsState << " ";

    return os;

}
std::istream& Spin::Read(std::istream& is)
{
    int temp;
    if (StreamableObject::Binary())
        BinaryRead(temp,is);
    else
        is >> temp;
    assert(is);
    itsState=static_cast<State>(temp);
    return is;
}

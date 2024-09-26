// File: Spin.C  A 3 state variable indicating the spin polarization.

#include "Misc/Spin.H"
#include "oml/imp/stream.h"
#include "oml/imp/binio.h"

std::ostream& operator<<(std::ostream& os,const Spin& S)
{
    if (StreamableObject::Binary())
    {
        int temp=S.itsState;
        BinaryWrite(temp,os);
    }
    else
        os << S.itsState << " ";

    return os;

}
std::istream& operator>>(std::istream& is, Spin& S)
{
    int temp;
    if (StreamableObject::Binary())
        BinaryRead(temp,is);
    else
        is >> temp;

    S.itsState=(Spin::State)temp;
    return is;
}

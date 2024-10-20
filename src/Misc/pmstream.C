#include "Imp/Misc/pmstream.h"
#include <typeinfo>

PMStreamableObject::~PMStreamableObject()
{

}

std::ostream& operator<<(std::ostream& os, const PMStreamableObject& o)
{
    o.WriteHeader(os,typeid(o).name());
    return o.Write(os);
}

std::istream& operator>>(std::istream& is,       PMStreamableObject& o)
{
    StreamableObject::Mode current=o.ReadHeader(is,typeid(o).name());
    o.Read (is);
    o.SetOutputMode(current); //Restore to previous state.
    return is;
}





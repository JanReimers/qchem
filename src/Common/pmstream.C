#include "Common/pmstream.h"
#include <typeinfo>

PMStreamableObject::~PMStreamableObject()
{

}

std::ostream& operator<<(std::ostream& os, const PMStreamableObject& o)
{
    // o.WriteHeader(os,typeid(o).name());
    return o.Write(os);
}


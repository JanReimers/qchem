#include "Base/stl_io.h"

//
//  W/R headers with output mode and typename
//
void WriteHeader(std::ostream& os,c_str type)
{
  assert(os);
  assert(type);
  if (!StreamableObject::Pretty())  os << (int)StreamableObject::GetOutputMode() << " " << type << " ";
  assert(type);
  assert(os);
}

StreamableObject::Mode ReadHeader(std::istream& is,c_str type)
{
    assert(is);
    assert(type);
    StreamableObject::Mode current,temp;
    int itemp;
    is >> itemp;
    temp=static_cast<StreamableObject::Mode>(itemp);
    assert(temp!=StreamableObject::pretty);
    current=StreamableObject::SetOutputMode(temp);
    StreamableObject::CheckName(is,type);
    assert(is);
    return current;
}

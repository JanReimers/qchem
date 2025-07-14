#ifndef PMSTREAM_H
#define PMSTREAM_H

#include <iosfwd>
//
// Provide op<< for any class that implements Write
//
class PMStreamableObject 
{
public: 
    virtual ~PMStreamableObject() {};
    virtual std::ostream& Write(std::ostream&) const=0;
};

inline std::ostream& operator<<(std::ostream& os, const PMStreamableObject& o)
{
    return o.Write(os);
}


#endif // PMSTREAM_H

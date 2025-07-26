// Common/Stream.C  Make classes streamable with op<<
module;
#include <iosfwd>

export module qchem.Streamable;
//
// Provide op<< for any class that implements Write
//
export class Streamable 
{
public: 
    virtual ~Streamable() {};
    virtual std::ostream& Write(std::ostream&) const=0;
};

export inline std::ostream& operator<<(std::ostream& os, const Streamable& o)
{
    return o.Write(os);
}
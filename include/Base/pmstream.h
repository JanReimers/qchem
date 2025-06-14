#ifndef PMSTREAM_H
#define PMSTREAM_H

#include "oml/imp/stream.h"
#include <iosfwd>

//
// Provide IO (Pickling) with virtual dispatch.
//
class PMStreamableObject : public StreamableObject
{
    public: 
        virtual ~PMStreamableObject();

        virtual std::ostream& Write(std::ostream&) const=0;
        virtual std::istream& Read (std::istream& is) {return is;}
};

std::ostream& operator<<(std::ostream& os, const PMStreamableObject& o);
std::istream& operator>>(std::istream& is,       PMStreamableObject& o);

#endif // PMSTREAM_H

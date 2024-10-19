// File: UniqueID.C  Anything derived from this will have a unique ID.



#include "Imp/Misc/UniqueID/UniqueID.H"
#include "oml/imp/binio.h"
#include <fstream>
#include <cassert>
#include <stdlib.h>

char ID_FILE[] = "/tmp/NextID.tmp";

UniqueID::IDtype UniqueID::NextID = 0      ;
const UniqueID::IDtype UniqueID::MaxID  = 0x20000;

UniqueID::UniqueID()
{
    if (!NextID)
    {
        if(!std::ifstream(ID_FILE))
        {
            std::ofstream Next(ID_FILE);
            Next << 1UL << std::endl;
        }
        std::ifstream Next(ID_FILE);
        if(Next)
            Next >> NextID;
        else
        {
            std::cerr << "Can't open NextID file" << std::endl;
            exit(-1);
        }

    }
    itsID=GetNextID();
}

UniqueID::IDtype UniqueID::GetNextID()
{
    if (++NextID >= GetMaxID()) NextID=1;
    return NextID;
}

UniqueID::UniqueID(const UniqueID&)
    : itsID(GetNextID())
{};

UniqueID::~UniqueID()
{
    std::ofstream Next(ID_FILE);
    Next << NextID << std::endl;
}

UniqueID& UniqueID::operator=(const UniqueID&)
{
    itsID=GetNextID();
    return *this;
}



std::ostream& UniqueID::Write(std::ostream& os) const
{
    if (StreamableObject::Binary()) BinaryWrite(itsID,os);
    if (StreamableObject::Ascii ()) os << itsID << " ";
    return os;
}

std::istream& UniqueID::Read (std::istream& is)
{
    if (StreamableObject::Binary()) BinaryRead(itsID,is);
    if (StreamableObject::Ascii ()) is >> itsID;
    return is;
}



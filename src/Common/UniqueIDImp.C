// File: UniqueIDImp.C  Anything derived from this will have a unique ID.



#include "Common/UniqueIDImp.H"
#include "oml/imp/binio.h"
#include <fstream>
#include <cassert>
#include <stdlib.h>

char ID_FILE[] = "/tmp/NextID.tmp";

UniqueID::IDtype UniqueIDImp::NextID = 0      ;
const UniqueID::IDtype UniqueIDImp::MaxID  = 0x20000;

UniqueIDImp::UniqueIDImp()
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

UniqueID::IDtype UniqueIDImp::GetNextID()
{
    if (++NextID >= GetMaxID()) NextID=1;
    return NextID;
}

UniqueIDImp::UniqueIDImp(const UniqueIDImp&)
    : itsID(GetNextID())
{};

UniqueIDImp::~UniqueIDImp()
{
    std::ofstream Next(ID_FILE);
    Next << NextID << std::endl;
}

UniqueID& UniqueIDImp::operator=(const UniqueID&)
{
    itsID=GetNextID();
    return *this;
}



std::ostream& UniqueIDImp::Write(std::ostream& os) const
{
    return os << itsID << " ";
}




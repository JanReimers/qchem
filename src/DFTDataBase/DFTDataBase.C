// File: DFTDataBase.C  Implementation for a calculation manager.



#include "DFTDataBase/DFTDataBase.H"
#include <iostream>
#include <fstream>
#include <stdlib.h>


DFTDataBase::DFTDataBase(const char* theFileName)
    : itsFileName(theFileName)
    , itsIntegralDataBases()
{
    if (itsFileName=="") itsFileName="DFTDataBase.tmp"; //Default file name
    std::ifstream in(itsFileName.c_str());
    if (!in)
    {
        std::cerr << "DFTDataBase::DFTDataBase Cannot open database file: " << theFileName << std::endl;
//        exit(-1);
    }
    else
    {
        //in >> itsIntegralDataBases;
        assert(in);
//        StreamableObject::Mode m=StreamableObject::SetToAscii();
//        std::ofstream fs("AtomAscii.db");
//        fs << itsIntegralDataBases;
//        fs.close();
//        StreamableObject::SetOutputMode(m);
    }
}

DFTDataBase::~DFTDataBase()
{
    StreamableObject::Mode mode=StreamableObject::SetToBinary();
    std::ofstream out(itsFileName.c_str());
    //out << itsIntegralDataBases;
    StreamableObject::SetOutputMode(mode);
}


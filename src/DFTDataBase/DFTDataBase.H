// File: DFTDataBase.H  Persistent data base.
#ifndef _DFTDataBase_H_
#define _DFTDataBase_H_



#include "Misc/Persistent/PerDB.H"
#include "IntegralDataBase.H"
#include <string>

//------------------------------------------------------------------------
//
//  This is a persistent data base for the OODFT system.  Upon destruction
//  the database is streamed out to disk.
//
class DFTDataBase
{
public:
    DFTDataBase(const char* FileName);
    ~DFTDataBase();

private:
    std::string                      itsFileName;
    PerDB<IntegralDataBase<double> > itsIntegralDataBases;
};


#endif

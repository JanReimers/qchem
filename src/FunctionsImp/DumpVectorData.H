// File: DumpVectorData.H  Dump plot data to a temp file.
#ifndef _DumpVectorData_H_
#define _DumpVectorData_H_



#include "Functions/VectorFunction.H"

class DumpVectorData :
    public virtual PlottableVectorFunction
{
public:
    DumpVectorData() {};
    virtual const char* DumpPlotData(const Mesh&, const RVec3& direction) const;
};

#endif //_DumpVectorData_H_

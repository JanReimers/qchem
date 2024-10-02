// File: DumpScalarData.H  Dump plot data to a temp file.
#ifndef _DumpScalarData_H_
#define _DumpScalarData_H_



#include "Functions/ScalarFunction.H"

class DumpScalarData :
    public virtual PlottableScalarFunction
{
public:
    DumpScalarData() {};
    virtual const char* DumpPlotData  (const Mesh&, const RVec3& direction  ) const;
    virtual const char* Dump3DPlotData(const Mesh&, const RVec3& n, double z) const;
};

#endif //_DumpScalarData_H_

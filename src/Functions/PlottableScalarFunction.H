// File: PlottableScalarFunction.H  Mixin interface for real space functions.
#ifndef _PlottableScalarFunction_H_
#define _PlottableScalarFunction_H_



#include "Functions/ScalarFunction.H"

class Mesh;

class PlottableScalarFunction
    : public virtual ScalarFunction<double>
{
public:
    virtual const char* DumpPlotData  (const Mesh&, const RVec3& direction  ) const=0; //Returns temporary file name.
    virtual const char* Dump3DPlotData(const Mesh&, const RVec3& n, double z) const=0; //Returns temporary file name.
};

#endif // _PlottableScalarFunction_H_

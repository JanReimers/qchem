// File: PlottableVectorFunction.H  Add plotting capability to vector functions.
#ifndef _PlottableVectorFunction_H_
#define _PlottableVectorFunction_H_



#include "Functions/VectorFunction.H"

class Mesh;

class PlottableVectorFunction
    : public virtual VectorFunction<double>
{
public:
    virtual const char* DumpPlotData (const Mesh&, const RVec3& direction) const=0; //Returns temporary data file name.
};

#endif // _PlottableVectorFunction_H_

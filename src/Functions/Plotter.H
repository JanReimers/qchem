// File: Plotter.H  Object the does the actual plotting.
#ifndef _Plotter_H_
#define _Plotter_H_

#include "types.H"


class Plotter
{
public:
    virtual void Plot   (const PlottableScalarFunction&, const Mesh*, const RVec3& direction)=0;
    virtual void Plot   (const PlottableVectorFunction&, const Mesh*, const RVec3& direction, int n=5)=0;
    virtual void AddPlot(const PlottableScalarFunction&)=0;
    virtual void Plot3D (const PlottableScalarFunction&, const Mesh*, const RVec3& n, double z=0)=0;

    virtual void LogXAxis()=0;
    virtual void LogYAxis()=0;
    virtual void LogZAxis()=0;
    virtual void Lines   ()=0;
    virtual void Symbols ()=0;
    virtual void Contours()=0;
};

#endif //_Plotter_H_

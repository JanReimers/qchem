// File: PlotterImplementation.H  Object the does the actual plotting.
#ifndef _PlotterImplementation_H_
#define _PlotterImplementation_H_



#include "Functions/Plotter.H"

#include "oml/vector3d.h"
#include <string>
#include <Misc/pstream.h>
#define PFSTREAM redi::opstream

typedef Vector3D<double> RVec;

//----------------------------------------------------------------------
//
//  Pipe all data and commands to Gnuplot.  Easy pie!
//
class PlotterImplementation
    : public virtual Plotter
{
public:
    PlotterImplementation();
    ~PlotterImplementation();

    virtual void Plot   (const PlottableScalarFunction&, const Mesh*, const RVec& direction);
    virtual void Plot   (const PlottableVectorFunction&, const Mesh*, const RVec& direction, int n=5);
    virtual void AddPlot(const PlottableScalarFunction&);
    virtual void Plot3D (const PlottableScalarFunction&, const Mesh*, const RVec& n, double z=0);

    virtual void LogXAxis();
    virtual void LogYAxis();
    virtual void LogZAxis();
    virtual void Lines   ();
    virtual void Symbols ();
    virtual void Contours();

private:
    PFSTREAM          ThePlotter;
    std::string       CurrentLineStyle;
    std::vector<std::string> itsDataFiles;
    const Mesh*       theMesh;
    RVec              theDirection;
};

#endif //_PlotterImplementation_H_

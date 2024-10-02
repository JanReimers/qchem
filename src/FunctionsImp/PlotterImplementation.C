// File: PlotterImplementation.C  Implementation for a plotter.



#include "FunctionsImp/PlotterImplementation.H"
#include "Functions/ScalarFunction.H"
#include "Functions/VectorFunction.H"
#include <unistd.h>
#include <cassert>

PlotterImplementation::PlotterImplementation()
    : ThePlotter("|gnuplot")
    , CurrentLineStyle("lines")
    , theMesh(0)
    , theDirection(1,0,0)
{
    ThePlotter << "load '.dftplot'" << std::endl;
};

PlotterImplementation::~PlotterImplementation()
{
    for (std::vector<std::string>::iterator i(itsDataFiles.begin()); i!=itsDataFiles.end(); i++)
        unlink(i->c_str());
    ThePlotter << "quit " << std::endl << std::endl;
}

void PlotterImplementation::Plot(const PlottableScalarFunction& func, const Mesh* X, const RVec& direction)
{
    theMesh=X;
    theDirection=direction;
    itsDataFiles.push_back(func.DumpPlotData(*theMesh,theDirection));
    ThePlotter << "plot '"<< itsDataFiles.back() << "' with " << CurrentLineStyle << std::endl;
}

void PlotterImplementation::AddPlot(const PlottableScalarFunction& func)
{
    assert(theMesh);
    itsDataFiles.push_back(func.DumpPlotData(*theMesh,theDirection));
    ThePlotter << "replot '"<< itsDataFiles.back() << "' with " << CurrentLineStyle << std::endl;
}

void PlotterImplementation::Plot(const PlottableVectorFunction& func, const Mesh* X, const RVec& direction, int n)
{
    theMesh=X;
    theDirection=direction;
    itsDataFiles.push_back(func.DumpPlotData(*theMesh,theDirection));
    ThePlotter << "plot ";
    for (index_t i=0; i<n-1; i++)
        ThePlotter << "'" << itsDataFiles.back() << "' using 1:" << i+2 << " with "  << CurrentLineStyle << ",";

    ThePlotter << "'" << itsDataFiles.back() << "' using 1:" << n+1 << " with "  << CurrentLineStyle << std::endl;
}

void PlotterImplementation::Plot3D(const PlottableScalarFunction& func, const Mesh* X, const RVec& n, double z)
{
    theMesh=X;
    theDirection=n;
    itsDataFiles.push_back(func.Dump3DPlotData(*theMesh,theDirection,z));
    ThePlotter << "splot '"<< itsDataFiles.back() << "' with lines" << std::endl;
}


void PlotterImplementation::LogXAxis()
{
    ThePlotter << "set logscale x" << std::endl;
}

void PlotterImplementation::LogYAxis()
{
    ThePlotter << "set logscale y" << std::endl;
}

void PlotterImplementation::LogZAxis()
{
    ThePlotter << "set logscale z" << std::endl;
}

void PlotterImplementation::Lines()
{
    CurrentLineStyle="lines";
}

void PlotterImplementation::Symbols()
{
    CurrentLineStyle="points";
}

void PlotterImplementation::Contours()
{
    ThePlotter << "set nosurface" << std::endl;
    ThePlotter << "set contour" << std::endl;
    ThePlotter << "set cntrparam levels auto 15" << std::endl;
    ThePlotter << "set view 0,90,1.5" << std::endl;
}


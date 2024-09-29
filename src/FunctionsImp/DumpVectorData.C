// File: DumpVectorData.C  Implementation for plotting.



#include "FunctionsImp/DumpVectorData.H"
#include "Mesh/Mesh.H"
#include "Mesh/MeshBrowser.H"
#include "oml/vector.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>

const char* DumpVectorData::DumpPlotData(const Mesh& X, const RVec3& direction) const
{
    RVec3 nd=normalize(direction);
    static char fname[256];
    char* ret=tmpnam(fname);
    (void)*ret; //Avoid unused warning.
    std::ofstream data(fname);
    MeshBrowser mb(X);
    for (; mb; mb++)
    {
        data << mb.R()*nd ;
        Vec y=(*this)(mb.R());
        Vec::const_iterator yb(y.begin());
        for (; yb!=y.end(); yb++) data << " " << *yb;
        data << std::endl;
    }
    data.close();
    return fname;
}


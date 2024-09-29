// File: DumpScalarData.cpp  Implementation for plotting.



#include "FunctionsImp/DumpScalarData.H"
#include "Mesh/Mesh.H"
#include "Mesh/MeshBrowser.H"
#include "oml/vector.h"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <unistd.h>
#include <cassert>
#include <cmath>


const char* DumpScalarData::DumpPlotData(const Mesh& X, const RVec3& direction) const
{
    RVec3 nd=normalize(direction);
    static char fname[256];
    char* ret=tmpnam(fname);
    (void)*ret;
    std::ofstream data(fname);
    MeshBrowser mb(X);
    for (; mb; mb++)
    {
        data << mb.R()*nd << " " << (*this)(mb.R()) << std::endl;
    }
    data.close();
    return fname;
}


void GetNormalCoordinates(const RVec3& z, RVec3& x, RVec3& y);

const char* DumpScalarData::Dump3DPlotData(const Mesh& X, const RVec3& n, double z) const
{
    assert(norm(n)>0);

    Vec  V(X.GetNumPoints());
    MeshBrowser m(X);
    for (index_t i=1; m; i++,m++) V(i)=m.R().x;
    RVec3 xhat,yhat;
    GetNormalCoordinates(n,xhat,yhat);

    static char fname[256];
    char* ret=tmpnam(fname);
    (void)*ret; //Avoid unused warning.
    std::ofstream data(fname);

    for (Vec::const_iterator x(V.begin()); x!=V.end(); x++)
    {
        for (Vec::const_iterator y(V.begin()); y!=V.end(); y++)
        {
            RVec3 r=*x*xhat + *y*yhat + z*n;
            data << (*this)(r) << std::endl;
        }
        data << std::endl;
    }
    data.close();
    return fname;
}

//
//  Given z find x and y such that x,y,z are orthnormal.
//
void GetNormalCoordinates(const RVec3& z, RVec3& xhat, RVec3& yhat)
{
    assert(norm(z)>0);
    RVec3 zhat=normalize(z);
//
// Deal with trivial cases where at least one component of z is 0.
// Some of these also represent singular cases which must be avoided
// for the general algoritm below.
//
    if (zhat.x==0)
    {
        xhat=RVec3(1,0,0);
        yhat=Cross(zhat,xhat);
        return;
    }
    if (zhat.y==0)
    {
        yhat=RVec3(0,1,0);
        xhat=Cross(yhat,zhat);
        return;
    }
    if (zhat.z==0)
    {
        xhat=RVec3(0,0,1);
        yhat=Cross(zhat,xhat);
        return;
    }
//
//  If we get this far then all components of z are non-zero.
//
    double d=sqrt(zhat.y*zhat.y + zhat.z*zhat.z);
    xhat=RVec3(0, -zhat.z/d, zhat.y/d);
    yhat=Cross(zhat,xhat);
    return;
};

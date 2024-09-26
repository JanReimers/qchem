// File: AngularMeshBrowser.H  Radial AngularMesh browser
#ifndef _AngularMeshBrowser_H_
#define _AngularMeshBrowser_H_



#include "oml/vector3d.h"
#include "oml/vector.h"
#include "oml/vector.h"

class AngularMesh;
typedef Vector3D<double> RVec3;

class AngularMeshBrowser
{
public:
    AngularMeshBrowser(const AngularMesh&);
    void operator++(int)
    {
        Db++;
        Wb++;
    }
    operator bool()
    {
        return Db!=DbEnd&&Wb!=WbEnd;
    }
    const RVec3  D() const
    {
        return *Db;
    }
    const double W() const
    {
        return *Wb;
    }

private:
    Vector<RVec3> ::const_iterator Db;
    Vector<double>::const_iterator Wb;
    Vector<RVec3> ::const_iterator DbEnd;
    Vector<double>::const_iterator WbEnd;
};

#endif //_AngularMeshBrowser_H_


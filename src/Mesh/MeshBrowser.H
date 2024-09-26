// File: MeshBrowser.H  Radial mesh browser
#ifndef _MeshBrowser_H_
#define _MeshBrowser_H_



#include "oml/vector3d.h"
#include "oml/vector.h"
#include "oml/vector.h"

class Mesh;
typedef Vector3D<double> RVec3;

class MeshBrowser
{
public:
    MeshBrowser(const Mesh&);
    void operator++(int)
    {
        Rb++;
        Wb++;
    }
    operator bool()
    {
        return Rb!=RbEnd&&Wb!=WbEnd;
    }
    const RVec3  R() const
    {
        return *Rb;
    }
    const double W() const
    {
        return *Wb;
    }

private:
    Vector<RVec3> ::const_iterator Rb;
    Vector<double>::const_iterator Wb;
    Vector<RVec3> ::const_iterator RbEnd;
    Vector<double>::const_iterator WbEnd;
};

#endif //_MeshBrowser_H_


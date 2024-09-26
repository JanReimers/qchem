// File: RadialMeshBrowser.H  Radial mesh browser
#ifndef _RadialMeshBrowser_H_
#define _RadialMeshBrowser_H_



#include "oml/vector.h"

class RadialMesh;

class RadialMeshBrowser
{
public:
    RadialMeshBrowser(const RadialMesh&);
    void operator++(int)
    {
        Rb++;
        Wb++;
    }
    operator bool()
    {
        return Rb!=RbEnd&&Wb!=WbEnd;
    }
    double R() const
    {
        return *Rb;
    }
    double W() const
    {
        return *Wb;
    }

private:
    Vector<double>::const_iterator Rb;
    Vector<double>::const_iterator Wb;
    Vector<double>::const_iterator RbEnd;
    Vector<double>::const_iterator WbEnd;
};

#endif //_RadialMeshBrowser_H_


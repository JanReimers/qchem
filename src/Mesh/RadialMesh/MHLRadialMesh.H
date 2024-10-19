// File: MHLRadialMesh.H  Murry, Handy Laming radial mesh.
#ifndef _MHLRadialMesh_H_
#define _MHLRadialMesh_H_



#include "Mesh/RadialMesh/RadialMesh.H"

//
// C.W. Murray,N.C. HandyandG.J. Laming, Mol. Phys.78 (1993) 997. 
//
class MHLRadialMesh : public  RadialMesh
{
public:
    MHLRadialMesh(index_t NumPoints, int m, double alpha);
};

#endif //_MHLRadialMesh_H_

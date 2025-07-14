// File: UnitCellMesh.H  UnitCell mesh, mostly used for plotting.
#ifndef _UnitCellMesh_H_
#define _UnitCellMesh_H_

import Mesh;

class UnitCell;

class UnitCellMesh : public Mesh
{
public:
    UnitCellMesh() {};
    UnitCellMesh(const UnitCell&, index_t numPoints);
    Mesh*  Clone() const;
};

#endif //_UnitCellRadialMesh_H_

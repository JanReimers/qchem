// File: UnitCellMesh.C  UnitCell mesh implementation.
module;
#include <cmath>

export module Cluster.UnitCellMesh;
import Cluster.UnitCell;
import qchem.Mesh;
import Cluster.UnitCell;
import oml;

export class UnitCellMesh : public Mesh
{
public:
    UnitCellMesh() {};
    UnitCellMesh(const UnitCell&, size_t numPoints);
    Mesh*  Clone() const;
};

UnitCellMesh::UnitCellMesh(const UnitCell& cell, size_t NumPoints)
{
    size_t N=NumPoints*NumPoints*NumPoints;
    double w=cell.GetCellVolume()/N;
    double del=1.0/NumPoints;

    for(size_t i1=0; i1<NumPoints; i1++)
        for(size_t i2=0; i2<NumPoints; i2++)
            for(size_t i3=0; i3<NumPoints; i3++)
                push_back(RVec3(i1*del,i2*del,i3*del),w);
}

Mesh* UnitCellMesh::Clone() const
{
    return new UnitCellMesh(*this);
}


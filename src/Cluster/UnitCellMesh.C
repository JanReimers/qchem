// File: UnitCellMesh.C  UnitCell mesh implementation.



#include "Cluster/UnitCell.H"
#include "UnitCellMesh.H"
#include <Common/Constants.H>
#include <cmath>

UnitCellMesh::UnitCellMesh(const UnitCell& cell, index_t NumPoints)
{
    index_t N=NumPoints*NumPoints*NumPoints;
    double w=cell.GetCellVolume()/N;
    double del=1.0/NumPoints;

    for(int i1=0; i1<NumPoints; i1++)
        for(int i2=0; i2<NumPoints; i2++)
            for(int i3=0; i3<NumPoints; i3++)
                push_back(RVec3(i1*del,i2*del,i3*del),w);
}

Mesh* UnitCellMesh::Clone() const
{
    return new UnitCellMesh(*this);
}


// File: UnitCellMesh.C  UnitCell mesh implementation.



#include "Imp/Cluster/UnitCell.H"
#include "Mesh/UnitCellMesh.H"
#include "Misc/DFTDefines.H"
#include <cmath>

UnitCellMesh::UnitCellMesh(const UnitCell& cell, index_t NumPoints)
    : MeshImplementation(NumPoints*NumPoints*NumPoints)
{
    index_t N=NumPoints*NumPoints*NumPoints;

    Vector<RVec3>  Rv(N);
    Vector<double> W (N);

    Fill(W,cell.GetCellVolume()/N);

    double del=1.0/NumPoints;
    Vector<RVec3>::iterator ri(Rv.begin());

    for(int i1=0; i1<NumPoints; i1++)
        for(int i2=0; i2<NumPoints; i2++)
            for(int i3=0; i3<NumPoints; i3++,ri++)
                *ri=RVec3(i1*del,i2*del,i3*del);

    Initialize(Rv,W);
}

Mesh* UnitCellMesh::Clone() const
{
    return new UnitCellMesh(*this);
}


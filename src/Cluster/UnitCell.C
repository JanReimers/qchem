// File: UnitCell.C  Unit cell for a lattice.
module;
#include "Common/pmstream.h"

export module Cluster.UnitCell;
import oml;

export class UnitCell
    : public virtual PMStreamableObject
{
public:
     using RVec3=Vector3D<double>;
    UnitCell();
    UnitCell(double a, double b, double c, double alpha, double beta, double gamma);
    UnitCell(const Matrix3D<double> MetricTensor);

    UnitCell      MakeReciprocalCell() const;

    double        GetCellVolume     (                  ) const;
    double        GetMinimumCellEdge(                  ) const;
    double        GetDistance       (const RVec3& r    ) const;
    Vector3D<int> GetNumCells       (double MaxDistance) const;

    std::ostream&  Write(std::ostream&) const;

    static UnitCell* Factory(std::istream&);

private:
    double itsA,itsB,itsC;              //Angstroms.
    double itsAlpha, itsBeta, itsGamma; //Radians.

    Matrix3D<double> itsMetricTensor;
};



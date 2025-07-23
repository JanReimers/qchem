// File: UnitCell.C  Unit cell for a lattice.
module;
#include <iosfwd>
export module Cluster.UnitCell;
export import qchem.Types;
import oml.Matrix3D;
import qchem.Streamable;

export class UnitCell
    : public virtual Streamable
{
public:
    UnitCell();
    UnitCell(double a, double b, double c, double alpha, double beta, double gamma);
    UnitCell(const Matrix3D<double> MetricTensor);

    UnitCell      MakeReciprocalCell() const;

    double        GetCellVolume     (                  ) const;
    double        GetMinimumCellEdge(                  ) const;
    double        GetDistance       (const RVec3& r    ) const;
    Vector3D<int> GetNumCells       (double MaxDistance) const;

    std::ostream&  Write(std::ostream&) const;

private:
    double itsA,itsB,itsC;              //Angstroms.
    double itsAlpha, itsBeta, itsGamma; //Radians.

    Matrix3D<double> itsMetricTensor;
};



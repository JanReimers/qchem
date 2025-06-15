// File: UnitCell.H  Unit cell for a lattice.
#ifndef _UnitCell_H_
#define _UnitCell_H_

#include "types.H"
#include "oml/vector3d.h"
#include "oml/matrix3d.h"
#include "Common/pmstream.h"

class UnitCell
    : public virtual PMStreamableObject
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
    std::istream&  Read (std::istream&)      ;

    static UnitCell* Factory(std::istream&);

private:
    double itsA,itsB,itsC;              //Angstroms.
    double itsAlpha, itsBeta, itsGamma; //Radians.

    Matrix3D<double> itsMetricTensor;
};

#endif // _UnitCell_H_



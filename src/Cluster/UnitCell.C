// File: UnitCell.C  Unit cell for a lattice.
module;
#include <cmath>
#include <iostream>
#include <cassert>
#include "Common/pmstream.h"

export module Cluster.UnitCell;
import Common.Constants;
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



UnitCell::UnitCell()
    : itsA    (10.0)
    , itsB    (10.0)
    , itsC    (10.0)
    , itsAlpha(Rad90)
    , itsBeta (Rad90)
    , itsGamma(Rad90)
    , itsMetricTensor
    (
        itsA*itsA              , itsB*itsA*cos(itsGamma), itsC*itsA*cos(itsBeta ),
        itsA*itsB*cos(itsGamma), itsB*itsB              , itsC*itsB*cos(itsAlpha),
        itsA*itsC*cos(itsBeta ), itsB*itsC*cos(itsAlpha), itsC*itsC
    )
{
}

UnitCell::UnitCell(double a, double b, double c, double alpha, double beta, double gamma)
    : itsA    (a)
    , itsB    (b)
    , itsC    (c)
    , itsAlpha(Rad(alpha))
    , itsBeta (Rad(beta ))
    , itsGamma(Rad(gamma))
    , itsMetricTensor
    (
        itsA*itsA              , itsB*itsA*cos(itsGamma), itsC*itsA*cos(itsBeta ),
        itsA*itsB*cos(itsGamma), itsB*itsB              , itsC*itsB*cos(itsAlpha),
        itsA*itsC*cos(itsBeta ), itsB*itsC*cos(itsAlpha), itsC*itsC
    )
{};

UnitCell::UnitCell(const Matrix3D<double> m)
    : itsA    (sqrt(m(1,1)))
    , itsB    (sqrt(m(2,2)))
    , itsC    (sqrt(m(3,3)))
    , itsAlpha(acos(m(3,2)/(itsB*itsC)))
    , itsBeta (acos(m(3,1)/(itsA*itsC)))
    , itsGamma(acos(m(2,1)/(itsA*itsB)))
    , itsMetricTensor(m)
{}

UnitCell UnitCell::MakeReciprocalCell() const
{
    return UnitCell(2*Pi*Invert(itsMetricTensor));
}

double UnitCell::GetCellVolume() const
{
    return sqrt(Determinant(itsMetricTensor));
}

double UnitCell::GetMinimumCellEdge() const
{
    if (itsA < itsB)
    {
        return itsA<itsC ? itsA : itsC;
    }
    else
    {
        return itsB<itsC ? itsB : itsC;
    }
}

Vector3D<int> UnitCell::GetNumCells(double MaxDistance) const
{
    assert(MaxDistance>0);
    return Vector3D<int>((int)ceil(MaxDistance/itsA),(int)ceil(MaxDistance/itsB),(int)ceil(MaxDistance/itsC));
}

double UnitCell::GetDistance(const RVec3& r) const
{
    double ret=sqrt(r*itsMetricTensor*r);
    return ret;
}

std::ostream&  UnitCell::Write(std::ostream& os) const
{
    os << "(a,b,c)=(" << itsA << "," << itsB << "," << itsC << "),   "
    << "(alpha,beta,gamma)=(" << itsAlpha << "," << itsBeta << "," << itsGamma << ")";
    return os;
}





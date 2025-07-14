// File: UnitCell.C  Unit cell for a lattice.

#include "Cluster/UnitCell.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include <cmath>
#include <iostream>
#include <cassert>

import Common.Constants;

const double Rad90=M_PI/2.0;

inline double Rad(double d)
{
    return d/180.0*M_PI;
}


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
    if(Binary())
    {
        BinaryWrite(itsA    ,os);
        BinaryWrite(itsB    ,os);
        BinaryWrite(itsC    ,os);
        BinaryWrite(itsAlpha,os);
        BinaryWrite(itsBeta ,os);
        BinaryWrite(itsGamma,os);
    }
    if(Ascii()) os << itsA << " " << itsB << " " << itsC << " "
        << itsAlpha << " " << itsBeta << " " << itsGamma << " ";
    if (!StreamableObject::Pretty()) os << itsMetricTensor;
    if (StreamableObject::Pretty())
    {
        os << "(a,b,c)=(" << itsA << "," << itsB << "," << itsC << "),   "
        << "(alpha,beta,gamma)=(" << itsAlpha << "," << itsBeta << "," << itsGamma << ")";
    }
    return os;
}

std::istream&  UnitCell::Read (std::istream& is)
{
    if(Binary())
    {
        BinaryRead(itsA    ,is);
        BinaryRead(itsB    ,is);
        BinaryRead(itsC    ,is);
        BinaryRead(itsAlpha,is);
        BinaryRead(itsBeta ,is);
        BinaryRead(itsGamma,is);
    }
    if(Ascii()) is >> itsA >> itsB >> itsC >> itsAlpha >> itsBeta >> itsGamma;
    is >> itsMetricTensor;
    return is;
}




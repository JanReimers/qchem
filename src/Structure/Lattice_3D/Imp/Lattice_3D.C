// File: Structure/Lattice_3D/Imp/Lattice_3D.C Define a 3D infinite lattice.
module;
#include <iostream>
#include <cassert>
#include <algorithm> //sort
#include <vector>
#include <memory>    //make_shared (GetStructure)

module qchem.Lattice_3D;

import qchem.Streamable;
import qchem.Math;
import qchem.Blaze; //for op!= on blaze iterators (range-for over rvec3vec_t)

//--------------------------------------------------------------------------
//
//  Construction zone.
//
// Lattice_3D::Lattice(            )
//     : itsUnitCell (            )
//     , itsLimits   (0,0,0       )
//     , itsAtoms    (new Molecule)
//     , itsTolerence(0.0001      )
// {}

Lattice_3D::Lattice_3D(const UnitCell& cell, const Vector3D<int>& Limits)
    : itsUnitCell (cell        )
    , itsLimits   (Limits      )
    , itsTolerence(0.0001      )
{}


//---------------------------------------------------------
//
//  Structure stuff.
//

ReciprocalLattice Lattice_3D::Reciprocal() const
{
    return ReciprocalLattice(itsUnitCell.MakeReciprocalCell());
}

std::shared_ptr<const Structure> Lattice_3D::GetStructure() const
{
    return std::make_shared<UnitCell>(itsUnitCell); // deep-copies the atom basis (Cartesian a.u.)
}


//----------------------------------------------------------
//
//  Simple lattice questions.
//
size_t Lattice_3D::GetNumSites() const
{
    return GetNumBasisSites() * GetNumUnitCells();
}

size_t Lattice_3D::GetNumBasisSites() const
{
    return itsUnitCell.GetNumAtoms();
}

size_t Lattice_3D::GetNumUnitCells() const
{
    return itsLimits.x * itsLimits.y * itsLimits.z;
}

//------------------------------------------------------------------
//
//  Coordinate to site number translations.
//
size_t Lattice_3D::GetSiteNumber(const rvec3_t& r) const
{
    rvec3_t basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    size_t ib = Find(basis); //Find within tolerence.
    assert (ib<GetNumBasisSites());

    size_t sitenum=ib + GetNumBasisSites()*(cell.z + itsLimits.z*(cell.y + itsLimits.y*cell.x));
    assert(sitenum<GetNumSites());
    // std::cout << "GetCoordinate(sitenum)="  << GetCoordinate(sitenum) << std::endl;
    // std::cout << "rvec3_t(cell.x,cell.y,cell.z)="  << rvec3_t(cell.x,cell.y,cell.z) << std::endl;
    // std::cout << "basis="  << basis << std::endl;
    assert(itsUnitCell.GetDistance(GetCoordinate(sitenum)-rvec3_t(cell.x,cell.y,cell.z)-basis) < itsTolerence);
    return sitenum;
}

size_t Lattice_3D::GetBasisNumber(const rvec3_t& r) const
{
    rvec3_t basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    size_t ret=Find(basis);
    assert(ret<GetNumBasisSites());
    return ret;
}

size_t Lattice_3D::GetBasisNumber(size_t SiteNumber) const
{
    assert(SiteNumber<GetNumSites());
    return SiteNumber%GetNumBasisSites();
}

void Lattice_3D::SplitCoordinate(const rvec3_t& r, rvec3_t& basis, Vector3D<int>& cell) const
{
    cell.x=(int)floor(r.x);
    cell.y=(int)floor(r.y);
    cell.z=(int)floor(r.z);

    basis.x=r.x-cell.x;
    basis.y=r.y-cell.y;
    basis.z=r.z-cell.z;

    cell.x = cell.x%itsLimits.x;
    cell.y = cell.y%itsLimits.y;
    cell.z = cell.z%itsLimits.z;

    if (cell.x < 0) cell.x+=itsLimits.x;
    if (cell.y < 0) cell.y+=itsLimits.y;
    if (cell.z < 0) cell.z+=itsLimits.z;
}

Vector3D<int> Lattice_3D::GetCellCoord (const rvec3_t& r) const
{
    rvec3_t basis;
    Vector3D<int> ret;
    SplitCoordinate(r,basis,ret);
    return ret;
}


rvec3_t Lattice_3D::GetCoordinate(size_t SiteNumber) const
{
    assert(SiteNumber<GetNumSites());
    size_t ib=GetBasisNumber(SiteNumber);
    rvec3_t ret=GetBasisVector(ib);

    SiteNumber-=ib;
    SiteNumber/=GetNumBasisSites();
    int iz=SiteNumber%itsLimits.z;
    assert(iz>=0);
    assert(iz<itsLimits.z);

    SiteNumber-=iz;
    SiteNumber/=itsLimits.z;
    int iy=SiteNumber%itsLimits.y;
    assert(iy>=0);
    assert(iy<itsLimits.y);

    SiteNumber-=iy;
    SiteNumber/=itsLimits.y;
    int ix=SiteNumber;
    assert(ix>=0);
    assert(ix<itsLimits.x);

    ret+=rvec3_t(ix,iy,iz);
    return ret;
}

//----------------------------------------------------------
//
//  Advanced lattice questions.
//
rvec_t Lattice_3D::GetDistances(size_t NumShells) const
{
    double maxd=itsUnitCell.GetMinimumCellEdge()*NumShells; //Initial guess.
    std::vector<double> distances; //scratch buffer: built by push_back, then sorted

    rvec3vec_t super_cells=GetSuperCells(maxd);

    for (auto a1:itsUnitCell)
        for (auto a2:itsUnitCell)
            for (auto& c:super_cells)
            {
                double d=itsUnitCell.GetDistance(c + a2->itsR - a1->itsR);
                if(d>0 && d<=maxd && Find(d,distances)==distances.size()) distances.push_back(d);
            }

    std::sort(distances.begin(),distances.end());
    assert(distances.size()>=NumShells); //guess radius minEdge*NumShells should be ample
    return rvec_t(std::min(NumShells,distances.size()),distances.data());
}

rvec3vec_t Lattice_3D::GetBonds(size_t BasisNumber, double Distance) const
{
    assert(BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<rvec3_t> ret; //scratch buffer: size not known up front
    rvec3_t rb=GetBasisVector(BasisNumber);
    rvec3vec_t super_cells=GetSuperCells(Distance);

    for (auto a:itsUnitCell)
        for (auto& c:super_cells)
        {
            rvec3_t bond = a->itsR + c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (fabs(mbond-Distance) < itsTolerence) ret.push_back(bond);
        }
    return rvec3vec_t(ret.size(),ret.data());
}

rvec3vec_t Lattice_3D::GetBondsInSphere(size_t BasisNumber, double Distance) const
{
    assert(BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<rvec3_t> ret; //scratch buffer: size not known up front
    rvec3_t rb=GetBasisVector(BasisNumber);
    rvec3vec_t super_cells=GetSuperCells(Distance);

    for (auto a:itsUnitCell)
        for (auto& c:super_cells)
        {
            rvec3_t bond = a->itsR + c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (mbond<Distance+itsTolerence) ret.push_back(bond);
        }
    return rvec3vec_t(ret.size(),ret.data());
}

std::vector<ivec3_t>  Lattice_3D::GetCellsInSphere(double rmax) const
{
    return itsUnitCell.CellsInSphere(rmax); //direct lattice vectors R within rmax
}
//--------------------------------------------------------
//
//  Private unitilities.
//
size_t  Lattice_3D::Find(const rvec3_t& r) const //Search within the primary unit cell.
{
    size_t ret=GetNumBasisSites();
    size_t i=0;
    for (auto a:itsUnitCell)
    {
        if (itsUnitCell.GetDistance(r - a->itsR) < itsTolerence)
        {
            ret=i;
            break;
        }
        i++;
    } 
    return ret;
}

size_t  Lattice_3D::Find(double r,const std::vector<double>& lis) const
{
    size_t ret=lis.size();
    size_t i=0;
    for (std::vector<double>::const_iterator b(lis.begin()); b!=lis.end(); b++,i++) if (fabs(r-*b) < itsTolerence)
        {
            ret=i;
            break;
        }
    return ret;
}

rvec3vec_t Lattice_3D::GetSuperCells(double MaxDistance) const
{
    Vector3D<int> nc=itsUnitCell.GetNumCells(MaxDistance);
    rvec3vec_t ret((2*nc.x+1)*(2*nc.y+1)*(2*nc.z+1)); //full box, size known up front
    size_t k=0;
    for (int ix=-nc.x; ix<=nc.x; ix++)
        for (int iy=-nc.y; iy<=nc.y; iy++)
            for (int iz=-nc.z; iz<=nc.z; iz++)
                ret[k++]=ivec3_t(ix,iy,iz);
    return ret;
}

rvec3_t Lattice_3D::GetBasisVector(size_t BasisNumber) const
{
    assert(BasisNumber<GetNumBasisSites());
    rvec3_t ret;
    {
        for (auto b:itsUnitCell)
        {
            if (BasisNumber==0)
            {
                ret=b->itsR;
                break;
            }
            BasisNumber--;
        }
    }
    return ret;
}

using std::endl;
//------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& Lattice_3D::Write(std::ostream& os) const
{
    os << itsUnitCell << endl << itsLimits;

    return os;
}





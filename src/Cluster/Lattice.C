// File: Lattice.C



#include "Cluster/Lattice.H"
#include "Cluster/Molecule.H"
#include "Cluster/ClusterBrowser.H"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
//#include "oml/heapsort.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm> //sort

//--------------------------------------------------------------------------
//
//  Construction zone.
//
Lattice::Lattice(            )
    : itsUnitCell (            )
    , itsLimits   (0,0,0       )
    , itsAtoms    (new Molecule)
    , itsTolerence(0.0001      )
{}

Lattice::Lattice(const UnitCell& cell, const Vector3D<int>& Limits)
    : itsUnitCell (cell        )
    , itsLimits   (Limits      )
    , itsAtoms    (new Molecule)
    , itsTolerence(0.0001      )
{}

Lattice::Lattice(const UnitCell& cell, const Vector3D<int>& Limits,Cluster* Atoms)
    : itsUnitCell (cell  )
    , itsLimits   (Limits)
    , itsAtoms    (Atoms )
    , itsTolerence(0.0001)
{}

Lattice::~Lattice()
{
    delete itsAtoms;
}

//---------------------------------------------------------
//
//  Cluster stuff.
//
void Lattice::Insert(Atom* atom)
{
    itsAtoms->Insert(atom);
}

int Lattice::GetNumAtoms() const
{
    return itsAtoms->GetNumAtoms();
}

int Lattice::GetNuclearCharge() const
{
    return itsAtoms->GetNuclearCharge();
}

double Lattice::GetNetCharge() const
{
    return itsAtoms->GetNetCharge();
}

double Lattice::GetNumElectrons() const
{
    return itsAtoms->GetNumElectrons();
}

ChargeDensity* Lattice::GetChargeDensity() const
{
    return itsAtoms->GetChargeDensity();
}


//----------------------------------------------------------
//
//  Simple lattice questions.
//
int Lattice::GetNumSites() const
{
    return GetNumBasisSites() * GetNumUnitCells();
}

int Lattice::GetNumBasisSites() const
{
    return itsAtoms->GetNumAtoms();
}

int Lattice::GetNumUnitCells() const
{
    return itsLimits.x * itsLimits.y * itsLimits.z;
}

//------------------------------------------------------------------
//
//  Coordinate to site number translations.
//
int Lattice::GetSiteNumber(const RVec3& r) const
{
    RVec3 basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    int ib = Find(basis); //Find within tolerence.
    assert (ib<GetNumBasisSites());

    int sitenum=ib + GetNumBasisSites()*(cell.z + itsLimits.z*(cell.y + itsLimits.y*cell.x));
    assert(sitenum<GetNumSites());
    assert(sitenum>=0);
    assert(itsUnitCell.GetDistance(GetCoordinate(sitenum)-RVec3(cell.x,cell.y,cell.z)-basis) < itsTolerence);
    return sitenum;
}

int Lattice::GetBasisNumber(const RVec3& r) const
{
    RVec3 basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    int ret=Find(basis);
    assert(ret<GetNumBasisSites());
    return ret;
}

int Lattice::GetBasisNumber(int SiteNumber) const
{
    assert(SiteNumber>=0);
    assert(SiteNumber<GetNumSites());
    return SiteNumber%GetNumBasisSites();
}

void Lattice::SplitCoordinate(const RVec3& r, RVec3& basis, Vector3D<int>& cell) const
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

Vector3D<int> Lattice::GetCellCoord (const RVec3& r) const
{
    RVec3 basis;
    Vector3D<int> ret;
    SplitCoordinate(r,basis,ret);
    return ret;
}


RVec3 Lattice::GetCoordinate(int SiteNumber) const
{
    assert(SiteNumber>=0);
    assert(SiteNumber<GetNumSites());
    int ib=GetBasisNumber(SiteNumber);
    RVec3 ret=GetBasisVector(ib);

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

    ret+=RVec3(ix,iy,iz);
    return ret;
}

//----------------------------------------------------------
//
//  Advanced lattice questions.
//
std::vector<double> Lattice::GetDistances(int NumShells) const
{
    double maxd=itsUnitCell.GetMinimumCellEdge()*NumShells; //Initial guess.
    std::vector<double> distances;

    std::vector<RVec3> super_cells=GetSuperCells(maxd);

    for (ClusterBrowser b1(*itsAtoms); b1; b1++)
        for (ClusterBrowser b2(*itsAtoms); b2; b2++)
            for (std::vector<RVec3>::const_iterator c(super_cells.begin()); c!=super_cells.end(); c++)
            {
                double d=itsUnitCell.GetDistance(*c + b2->itsR - b1->itsR);
                if(d>0 && d<=maxd && Find(d,distances)==distances.size()) distances.push_back(d);
            }

    std::sort(distances.begin(),distances.end());
    return std::vector<double>(distances.begin(),distances.begin()+NumShells);
}

std::vector<RVec3> Lattice::GetBonds(int BasisNumber, double Distance) const
{
    assert(BasisNumber>=0 && BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<RVec3> ret;
    RVec3 rb=GetBasisVector(BasisNumber);
    std::vector<RVec3> super_cells=GetSuperCells(Distance);

    for (ClusterBrowser b(*itsAtoms); b; b++)
        for (std::vector<RVec3>::const_iterator c(super_cells.begin()); c!=super_cells.end(); c++)
        {
            RVec3 bond = b->itsR + *c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (fabs(mbond-Distance) < itsTolerence) ret.push_back(bond);
        }
    return ret;
}

std::vector<RVec3> Lattice::GetBondsInSphere(int BasisNumber, double Distance) const
{
    assert(BasisNumber>=0 && BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<RVec3> ret;
    RVec3 rb=GetBasisVector(BasisNumber);
    std::vector<RVec3> super_cells=GetSuperCells(Distance);

    for (ClusterBrowser b(*itsAtoms); b; b++)
        for (std::vector<RVec3>::const_iterator c(super_cells.begin()); c!=super_cells.end(); c++)
        {
            RVec3 bond = b->itsR + *c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (mbond<Distance+itsTolerence) ret.push_back(bond);
        }
    return ret;
}


//--------------------------------------------------------
//
//  Private unitilities.
//
index_t  Lattice::Find(const RVec3& r) const //Search within the primary unit cell.
{
    index_t ret=GetNumBasisSites();
    index_t i=0;
    for (ClusterBrowser b(*itsAtoms); b; b++,i++) if (itsUnitCell.GetDistance(r - b->itsR) < itsTolerence)
        {
            ret=i;
            break;
        }
    return ret;
}

index_t  Lattice::Find(double r,const std::vector<double>& lis) const
{
    index_t ret=lis.size();
    index_t i=0;
    for (std::vector<double>::const_iterator b(lis.begin()); b!=lis.end(); b++,i++) if (fabs(r-*b) < itsTolerence)
        {
            ret=i;
            break;
        }
    return ret;
}

std::vector<RVec3> Lattice::GetSuperCells(double MaxDistance) const
{
    std::vector<RVec3> ret;
    Vector3D<int> nc=itsUnitCell.GetNumCells(MaxDistance);
    for (int ix=-nc.x; ix<=nc.x; ix++)
        for (int iy=-nc.y; iy<=nc.y; iy++)
            for (int iz=-nc.z; iz<=nc.z; iz++)
                ret.push_back(RVec3(ix,iy,iz));
    return ret;
}

RVec3 Lattice::GetBasisVector(int BasisNumber) const
{
    assert(BasisNumber<GetNumBasisSites());
    RVec3 ret;
    {
        ClusterBrowser b(*itsAtoms);
        for (int ib=0; ib<BasisNumber; ib++,b++);
        ret=b->itsR;
    }
    return ret;
}

using std::endl;
//------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& Lattice::Write(std::ostream& os) const
{
    if ( Binary())
    {
        os << itsUnitCell << itsLimits << itsAtoms;
        BinaryWrite(itsTolerence,os);
    }
    if ( Ascii ()) os << itsUnitCell << itsLimits << endl << itsAtoms << endl << itsTolerence << " ";
    if ( Pretty()) os << itsUnitCell << endl << itsLimits << endl << itsAtoms << std::endl;
    return os;
}

std::istream& Lattice::Read (std::istream& is)
{
    is >> itsUnitCell >> itsLimits;
    delete itsAtoms;
    itsAtoms=Cluster::Factory(is);
    is >> itsAtoms;
    if (Binary())	BinaryRead(itsTolerence,is);
    if (Ascii ())	is >> itsTolerence;

    return is;
}



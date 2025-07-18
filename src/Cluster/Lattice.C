// File: Lattice.C
module;
#include <iostream>
#include <cassert>
#include <cmath>
#include <algorithm> //sort
#include <vector>
#include "MoleculeMesh.H"

export module Cluster.Lattice;

import Common.UniqueIDImp;
import qchem.Cluster;
import qchem.Atom;
import qchem.Molecule;
import Cluster.UnitCell;
import oml;

//----------------------------------------------------------------------------
//
//  A lattice is a big square block of unit cells.  Each unit cell has a bunch
//  of atoms in it.  You can create a molecule type cluster and insert it
//  into the lattice at construction time.  The molecule is then repated in each
//  unit cell.  Or you can create a lattice with no atoms.  Each atom inserted
//  goes into every unit cell.  For a bond structure calculation the spatial
//  extent of the lattice determines the number of k-points in the first zone,
//  and also determines the normalization constant for plane wave basis
//  functions.
//  ClustrerBrowsing through the lattice only really browses through the
//  primary unit cell.
//
export class Lattice
    : public virtual Cluster
    , public UniqueIDImp
{
public:
    Lattice();
    Lattice(const UnitCell&, const Vector3D<int>&);                //Empty unit cell.
    Lattice(const UnitCell&, const Vector3D<int>&,Cluster* Atoms); //Full  unit cell.
    ~Lattice();

    virtual void   Insert        (Atom*)      ;
    virtual size_t GetNumAtoms        () const;
    virtual int    GetNuclearCharge   () const;
    virtual double GetNetCharge       () const;
    virtual double GetNumElectrons    () const;
    virtual Mesh*  CreateMesh(const MeshParams&) const;

    //virtual ChargeDensity* GetChargeDensity   () const;

    const  UnitCell& GetUnitCell() const
    {
        return itsUnitCell;
    }

    double GetLatticeVolume() const
    {
        return GetNumSites()*itsUnitCell.GetCellVolume();
    }
    size_t    GetNumSites     () const;
    size_t    GetNumBasisSites() const;
    size_t    GetNumUnitCells () const;

    size_t    GetSiteNumber   (const RVec3&  ) const;
    size_t    GetBasisNumber  (const RVec3&  ) const;
    size_t    GetBasisNumber  (size_t SiteNumber) const;

    Vector3D<int> GetCellCoord   (const RVec3&  ) const;
    RVec3         GetCoordinate  (size_t SiteNumber) const;
    void          SplitCoordinate(const RVec3& r, RVec3& basis, Vector3D<int>& cell) const;

    std::vector<double> GetDistances    (size_t NumShells) const;
    std::vector<RVec3>  GetBonds        (size_t BasisNumber, double distance) const;
    std::vector<RVec3>  GetBondsInSphere(size_t BasisNumber, double distance) const;

    virtual const_iterator begin() const {return itsAtoms->begin();}
    virtual const_iterator end  () const {return itsAtoms->end  ();} 
    virtual       iterator begin()       {return itsAtoms->begin();}
    virtual       iterator end  ()       {return itsAtoms->end  ();} 
    
    std::ostream& Write(std::ostream&) const;

private:
    size_t      Find(const RVec3&               ) const; //Search unit cell.
    size_t      Find(double,const std::vector<double>&) const;
    std::vector<RVec3>  GetSuperCells(double MaxDistance) const;
    RVec3        GetBasisVector(size_t BasisNumber  ) const;

    UnitCell       itsUnitCell;  //Unit cell dimensions, no atoms.
    Vector3D<int>  itsLimits;    //Number of unit cell in each direction.
    Cluster*       itsAtoms;     //List of atoms in the unit cell.
    double         itsTolerence; //Positions closer than this are considered the same.
};

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

size_t Lattice::GetNumAtoms() const
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

Mesh*  Lattice::CreateMesh(const MeshParams& mp) const
{
    return new MoleculeMesh(*itsAtoms,mp);
}

//ChargeDensity* Lattice::GetChargeDensity() const
//{
//    return itsAtoms->GetChargeDensity();
//}


//----------------------------------------------------------
//
//  Simple lattice questions.
//
size_t Lattice::GetNumSites() const
{
    return GetNumBasisSites() * GetNumUnitCells();
}

size_t Lattice::GetNumBasisSites() const
{
    return itsAtoms->GetNumAtoms();
}

size_t Lattice::GetNumUnitCells() const
{
    return itsLimits.x * itsLimits.y * itsLimits.z;
}

//------------------------------------------------------------------
//
//  Coordinate to site number translations.
//
size_t Lattice::GetSiteNumber(const RVec3& r) const
{
    RVec3 basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    size_t ib = Find(basis); //Find within tolerence.
    assert (ib<GetNumBasisSites());

    size_t sitenum=ib + GetNumBasisSites()*(cell.z + itsLimits.z*(cell.y + itsLimits.y*cell.x));
    assert(sitenum<GetNumSites());
    assert(sitenum>=0);
    assert(itsUnitCell.GetDistance(GetCoordinate(sitenum)-RVec3(cell.x,cell.y,cell.z)-basis) < itsTolerence);
    return sitenum;
}

size_t Lattice::GetBasisNumber(const RVec3& r) const
{
    RVec3 basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    size_t ret=Find(basis);
    assert(ret<GetNumBasisSites());
    return ret;
}

size_t Lattice::GetBasisNumber(size_t SiteNumber) const
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


Lattice::RVec3 Lattice::GetCoordinate(size_t SiteNumber) const
{
    assert(SiteNumber>=0);
    assert(SiteNumber<GetNumSites());
    size_t ib=GetBasisNumber(SiteNumber);
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
std::vector<double> Lattice::GetDistances(size_t NumShells) const
{
    double maxd=itsUnitCell.GetMinimumCellEdge()*NumShells; //Initial guess.
    std::vector<double> distances;

    std::vector<RVec3> super_cells=GetSuperCells(maxd);

    for (auto& a1:*itsAtoms)
        for (auto& a2:*itsAtoms)
            for (std::vector<RVec3>::const_iterator c(super_cells.begin()); c!=super_cells.end(); c++)
            {
                double d=itsUnitCell.GetDistance(*c + a2->itsR - a1->itsR);
                if(d>0 && d<=maxd && Find(d,distances)==distances.size()) distances.push_back(d);
            }

    std::sort(distances.begin(),distances.end());
    return std::vector<double>(distances.begin(),distances.begin()+NumShells);
}

std::vector<Lattice::RVec3> Lattice::GetBonds(size_t BasisNumber, double Distance) const
{
    assert(BasisNumber>=0 && BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<RVec3> ret;
    RVec3 rb=GetBasisVector(BasisNumber);
    std::vector<RVec3> super_cells=GetSuperCells(Distance);

    for (auto& a:*itsAtoms)
        for (std::vector<RVec3>::const_iterator c(super_cells.begin()); c!=super_cells.end(); c++)
        {
            RVec3 bond = a->itsR + *c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (fabs(mbond-Distance) < itsTolerence) ret.push_back(bond);
        }
    return ret;
}

std::vector<Lattice::RVec3> Lattice::GetBondsInSphere(size_t BasisNumber, double Distance) const
{
    assert(BasisNumber>=0 && BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<RVec3> ret;
    RVec3 rb=GetBasisVector(BasisNumber);
    std::vector<RVec3> super_cells=GetSuperCells(Distance);

    for (auto& a:*itsAtoms)
        for (std::vector<RVec3>::const_iterator c(super_cells.begin()); c!=super_cells.end(); c++)
        {
            RVec3 bond = a->itsR + *c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (mbond<Distance+itsTolerence) ret.push_back(bond);
        }
    return ret;
}


//--------------------------------------------------------
//
//  Private unitilities.
//
size_t  Lattice::Find(const RVec3& r) const //Search within the primary unit cell.
{
    size_t ret=GetNumBasisSites();
    size_t i=0;
    for (auto& a:*itsAtoms)
    {
        i++;
        if (itsUnitCell.GetDistance(r - a->itsR) < itsTolerence)
        {
            ret=i;
            break;
        }
    } 
    return ret;
}

size_t  Lattice::Find(double r,const std::vector<double>& lis) const
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

std::vector<Lattice::RVec3> Lattice::GetSuperCells(double MaxDistance) const
{
    std::vector<RVec3> ret;
    Vector3D<int> nc=itsUnitCell.GetNumCells(MaxDistance);
    for (int ix=-nc.x; ix<=nc.x; ix++)
        for (int iy=-nc.y; iy<=nc.y; iy++)
            for (int iz=-nc.z; iz<=nc.z; iz++)
                ret.push_back(RVec3(ix,iy,iz));
    return ret;
}

Lattice::RVec3 Lattice::GetBasisVector(size_t BasisNumber) const
{
    assert(BasisNumber<GetNumBasisSites());
    RVec3 ret;
    {
        for (auto& b:*this)
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
std::ostream& Lattice::Write(std::ostream& os) const
{
    // os << itsUnitCell << endl << itsLimits << endl << itsAtoms << std::endl;
    return os;
}




// File: Lattice.h  Lattice of sites with basis.
#ifndef _Lattice_H_
#define _Lattice_H_

#include "Cluster.H"
#include "Imp/Cluster/UnitCell.H"
#include "Imp/Misc/UniqueID/UniqueIDImp.H"
#include "oml/vector3d.h"
#include <vector>

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
class Lattice
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
    std::istream& Read (std::istream&)      ;

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

#endif // _Lattice_H_

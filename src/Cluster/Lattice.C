// File: Lattice.C Define a 3D infite lattice.
module;
#include <vector>
#include <iosfwd>
#include <memory>

export module qchem.Lattice;
import Common.UniqueIDImp;
export import qchem.Cluster;
export import qchem.Atom;
export import Cluster.UnitCell;
import qchem.Mesh;

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
    // Lattice();
    typedef std::shared_ptr<Cluster> cl_t;
    Lattice(const UnitCell&, const Vector3D<int>&);                //Empty unit cell.
    Lattice(const UnitCell&, const Vector3D<int>&,const cl_t& Atoms); //Full  unit cell.

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
    Lattice   Reciprocal      (double Emax) const;  //Create the assosiated  reciprical Lattice;
    std::vector<rvec3_t>  GetReciprocalGrid() const;
    ivec3_t     GetLimits() const {return itsLimits;}

    size_t    GetNumSites     () const;
    size_t    GetNumBasisSites() const;
    size_t    GetNumUnitCells () const;

    size_t    GetSiteNumber   (const rvec3_t&  ) const;
    size_t    GetBasisNumber  (const rvec3_t&  ) const;
    size_t    GetBasisNumber  (size_t SiteNumber) const;

    Vector3D<int> GetCellCoord   (const rvec3_t&  ) const;
    rvec3_t         GetCoordinate  (size_t SiteNumber) const;
    void          SplitCoordinate(const rvec3_t& r, rvec3_t& basis, Vector3D<int>& cell) const;

    std::vector<double> GetDistances    (size_t NumShells) const;
    std::vector<rvec3_t>  GetBonds        (size_t BasisNumber, double distance) const;
    std::vector<rvec3_t>  GetBondsInSphere(size_t BasisNumber, double distance) const;
    std::vector<ivec3_t>  GetCellsInSphere(double distance);

    virtual const_iterator begin() const {return itsAtoms->begin();}
    virtual const_iterator end  () const {return itsAtoms->end  ();} 
    virtual       iterator begin()       {return itsAtoms->begin();}
    virtual       iterator end  ()       {return itsAtoms->end  ();} 
    
    std::ostream& Write(std::ostream&) const;

private:
    size_t      Find(const rvec3_t&               ) const; //Search unit cell.
    size_t      Find(double,const std::vector<double>&) const;
    std::vector<rvec3_t>  GetSuperCells(double MaxDistance) const;
    rvec3_t        GetBasisVector(size_t BasisNumber  ) const;

    UnitCell       itsUnitCell;  //Unit cell dimensions, no atoms.
    Vector3D<int>  itsLimits;    //Number of unit cell in each direction.
    cl_t           itsAtoms;     //List of atoms in the unit cell.
    double         itsTolerence; //Positions closer than this are considered the same.
};


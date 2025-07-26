// File: Cluster.C  A cluster of atoms, Molecule, unit cell, solid all with spatial struture defined.
module;
#include <vector>
#include <memory>

export module qchem.Cluster;
import Common.UniqueID; 
import qchem.Mesh; 
import qchem.Atom;
import oml.Vector3D;
import qchem.Streamable;

export class Cluster
    : public virtual UniqueID
    , public virtual Streamable
{
public:
    typedef std::vector<std::unique_ptr<Atom>> av_t;

    typedef av_t::      iterator       iterator;
    typedef av_t::const_iterator const_iterator;

    virtual ~Cluster() {};
    virtual bool operator==(const Cluster&) const
    {
        return false;
    }

    virtual void   Insert    (Atom*  )      =0;
    virtual size_t GetNumAtoms      () const=0;
    virtual int    GetNuclearCharge () const=0;
    virtual double GetNetCharge     () const=0;
    virtual double GetNumElectrons  () const=0;

    virtual size_t GetAtomIndex(const RVec3&, double tol=0.0) const;
    virtual Mesh*  CreateMesh(const MeshParams&) const=0;
    
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;

   
};





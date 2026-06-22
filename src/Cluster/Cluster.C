// File: Cluster.C  A structure of atoms, Molecule, unit cell, solid all with spatial struture defined.
module;
#include <vector>
#include <memory>

export module qchem.Cluster;
import Common.UniqueID; 
import qchem.Mesh; 
import qchem.Streamable;
import Common.UniqueIDImp;

export class Atom;

//----------------------------------------------------------------------------------------------
//! We an abstract interface that works for Atoms, Moleculas and lattice structures.
//! THis interface is used to implement Hamiltonian potentials for all three structures.
//! We need to be able to iterate through a list of atoms.  For a single atom structure
//! this means faking the iteration, which is a bit messy.
//! Each structure needs to know how to create an efficient numerical integration mesh.
//
export class Cluster
    : public virtual UniqueID
    , public virtual Streamable
{
public:
    typedef std::vector<Atom*> av_t;
    typedef av_t::      iterator       iterator;
    typedef av_t::const_iterator const_iterator;

    virtual ~Cluster() {};

    virtual std::string ID          () const;
    virtual size_t GetNumAtoms      () const=0;
    virtual int    GetNuclearCharge () const;
    virtual double GetNumElectrons  () const;
    virtual double GetNetCharge     () const;

    virtual size_t GetAtomIndex(const rvec3_t&, double tol=0.0) const;
    virtual Mesh*  CreateMesh(const MeshParams&) const=0;
     
    virtual const_iterator begin() const=0;
    virtual const_iterator end  () const=0;
    virtual       iterator begin()      =0;
    virtual       iterator end  ()      =0;
};

export class Atom
    : public Cluster
    , public UniqueIDImp
{
public:
    Atom(const Atom& a);
    Atom(int Z); 
    Atom(int Z, double charge);
    Atom(int Z, const rvec3_t& R);
    Atom(int Z, double charge, const rvec3_t& R);

    virtual size_t GetNumAtoms      () const {return 1;}
    virtual Mesh*  CreateMesh(const MeshParams&) const;

    virtual std::string   ID     () const;
    virtual std::ostream& Write  (std::ostream&) const;
    
    virtual const_iterator begin() const {return dummy.begin();}
    virtual const_iterator end  () const {return dummy.end  ();} 
    virtual       iterator begin()       {return dummy.begin();}
    virtual       iterator end  ()       {return dummy.end  ();} 

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    rvec3_t itsR;      //Spatial position.
private:
    av_t  dummy; //Kludge to get iterators from Atom* .
    Atom& operator=(const Atom&); //Why?  Rule of 5/6
};

export class Molecule 
    : public virtual Cluster
    , public UniqueIDImp
{
public:
    Molecule() {};
    Molecule(const Cluster& m);
    Molecule(const Molecule& m) : Molecule(static_cast<const Cluster&>(m)) {};
    virtual ~Molecule();
    virtual void   Insert     (Atom* a);
    virtual size_t GetNumAtoms() const;
    virtual Mesh*  CreateMesh (const MeshParams&) const;
       
    virtual std::ostream& Write(std::ostream&) const;

    virtual const_iterator begin() const {return itsAtoms.begin();}
    virtual const_iterator end  () const {return itsAtoms.end  ();} 
    virtual       iterator begin()       {return itsAtoms.begin();}
    virtual       iterator end  ()       {return itsAtoms.end  ();} 

private:
    av_t    itsAtoms;
    
};



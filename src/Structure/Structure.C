// File: Structure.C  A structure of atoms, Molecule, unit cell, solid all with spatial struture defined.
module;
#include <vector>
#include <memory>

export module qchem.Structure;
import qchem.Mesh;
import qchem.Streamable;
import Common.Iterators;

export class Atom;

//----------------------------------------------------------------------------------------------
//! We an abstract interface that works for Atoms, Moleculas and lattice structures.
//! THis interface is used to implement Hamiltonian potentials for all three structures.
//! We need to be able to iterate through a list of atoms.  For a single atom structure
//! this means faking the iteration, which is a bit messy.
//! Each structure needs to know how to create an efficient numerical integration mesh.
//
export class Structure
    : public virtual Streamable
{
public:
    virtual ~Structure() {};

    virtual std::string ID          () const;
    virtual size_t GetNumAtoms      () const=0;
    virtual int    GetNuclearCharge () const;
    virtual double GetNumElectrons  () const;
    virtual double GetNetCharge     () const;

    virtual size_t GetAtomIndex(const rvec3_t&, double tol=0.0) const;
    virtual Mesh*  CreateMesh(const MeshParams&) const=0;

    //  Element access and range-based iteration are built on the single
    //  storage primitive GetAtom() below, so a concrete Structure is free to
    //  store its atoms however it likes (Molecule keeps a vector, an Atom just
    //  returns itself, a future periodic cell could synthesize them).
    Atom* operator[](size_t i) const {return GetAtom(i);}
    typedef IndexIterator<Structure> const_iterator;
    const_iterator begin() const {return const_iterator(this,0);}
    const_iterator end  () const {return const_iterator(this,GetNumAtoms());}

protected:
    virtual Atom* GetAtom(size_t) const=0; //The only storage-specific primitive.
};

export class Atom
    : public Structure
{
public:
    Atom(int Z);
    Atom(int Z, double charge);
    Atom(int Z, const rvec3_t& R);
    Atom(int Z, double charge, const rvec3_t& R);
    //  Copy/move/assign are all defaulted: an Atom is now a plain value type
    //  (no self-reference, no UniqueID), so it is freely relocatable.

    virtual size_t GetNumAtoms      () const {return 1;}
    virtual Mesh*  CreateMesh(const MeshParams&) const;

    virtual std::string   ID     () const;
    virtual std::ostream& Write  (std::ostream&) const;

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    rvec3_t itsR;      //Spatial position.

protected:
    //  An Atom is a one-element structure consisting of itself.
    virtual Atom* GetAtom(size_t) const {return const_cast<Atom*>(this);}
};

export class Molecule
    : public virtual Structure
{
public:
    Molecule() {};
    Molecule(const Structure& m);
    Molecule(const Molecule& m) : Molecule(static_cast<const Structure&>(m)) {};
    virtual ~Molecule();
    virtual void   Insert     (Atom* a);
    virtual size_t GetNumAtoms() const;
    virtual Mesh*  CreateMesh (const MeshParams&) const;

    virtual std::ostream& Write(std::ostream&) const;

protected:
    virtual Atom* GetAtom(size_t i) const {return itsAtoms[i];}

private:
    std::vector<Atom*> itsAtoms;
};



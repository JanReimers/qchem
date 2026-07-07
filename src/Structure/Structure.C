// File: Structure.C  A structure of atoms, Molecule, unit cell, solid all with spatial struture defined.
module;
#include <vector>
#include <memory>
#include <functional>

export module qchem.Structure;
import qchem.Types;
import qchem.Streamable;
import qchem.Iterators;
import qchem.Mesh;        // qcMesh::Mesh / MeshParams -- the real-space integration mesh a geometry produces

namespace qchem {

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

    //! \brief True for a bounded (molecular/atomic) structure, false for a periodically-repeated one
    //! (a UnitCell/lattice).  The ion-ion (IonIon) Hamiltonian term uses this to choose a direct pair sum
    //! for a finite structure vs an Ewald lattice sum for a periodic one.
    virtual bool   isFinite         () const {return true;}

    //! \brief Sum of a per-element form factor over the structure's atoms, \f$\sum_a f(Z_a)\f$ -- a fair
    //! geometric question for ANY structure.  A periodic UnitCell normalises by its cell volume
    //! (\f$\frac1\Omega\sum_a f(Z_a)\f$, the G=0 background being a per-volume quantity), so a caller needs no
    //! cell-volume getter -- keeping \c CellVolume() off finite structures (no LSP violation).  The physics
    //! decision of WHETHER the sum contributes (e.g. the PP G=0 alignment is zero for a finite structure) is
    //! the caller's, via isFinite() -- this method only reports the geometric sum each structure can honestly give.
    virtual double SumFormFactors   (const std::function<double(int Z)>& f) const;

    //! \brief The real-space integration mesh natural to this geometry, at resolution \a mp.  The mesh TYPE
    //! is the STRUCTURE's responsibility, not the basis's: an Atom builds its single-centre (radial x angular)
    //! grid, a Molecule the multi-centre Becke grid, a periodic lattice a uniform / unit-cell-Becke grid.
    //! Pure virtual -- each geometry owns its most efficient mesh (no central dispatch).  A real-space
    //! orbital basis integrates a field on this via qcMesh::WeightedOverlap; plane waves own their G-grid and
    //! do not use it.
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const=0;

    virtual size_t GetAtomIndex(const rvec3_t&, double tol=0.0) const;

    //  Element access and range-based iteration are built on the single
    //  storage primitive GetAtom() below, so a concrete Structure is free to
    //  store its atoms however it likes (Molecule keeps a vector, an Atom just
    //  returns itself, a future periodic cell could synthesize them).
    Atom* operator[](size_t i) const {return GetAtom(i);}
    typedef IndexIterator<const Structure> const_iterator;
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

    virtual size_t GetNumAtoms      () const override {return 1;}

    //! A single atom's integration mesh: the atom-centred radial x angular grid (no Becke partitioning).
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    virtual std::string   ID     () const override;
    virtual std::ostream& Write  (std::ostream&) const override;

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    rvec3_t itsR;      //Spatial position.

protected:
    //  An Atom is a one-element structure consisting of itself.
    virtual Atom* GetAtom(size_t) const override {return const_cast<Atom*>(this);}
};

export class Molecule
    : public virtual Structure
{
public:
    Molecule() {};
    Molecule(const Structure& m);
    Molecule(const Molecule& m) : Molecule(static_cast<const Structure&>(m)) {};
    //  Owns its Atoms, so it needs a deep-copy assignment (rule of three);
    //  copy-and-swap reuses the copy ctor and is self-assignment safe.
    Molecule& operator=(Molecule m) {itsAtoms.swap(m.itsAtoms); return *this;}
    virtual ~Molecule();
    virtual void   Insert     (Atom* a);
    virtual size_t GetNumAtoms() const override;

    //! A molecule's integration mesh: the multi-centre Becke-partitioned grid.
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    virtual std::ostream& Write(std::ostream&) const override;

protected:
    virtual Atom* GetAtom(size_t i) const override {return itsAtoms[i];}

private:
    std::vector<Atom*> itsAtoms;
};

// Becke fuzzy-Voronoi integration mesh for a finite structure -- the shared implementation behind
// Atom::/Molecule::CreateIntegrationMesh.  Module-INTERNAL (NOT exported): callers use CreateIntegrationMesh,
// the public per-geometry entry point; this is the algorithm both geometries delegate to.  See Imp/Molecule.C
// (Becke 1988; coincident atoms R_ab=0 -> mu=0 so a coincident dimer integrates to the single-atom result).
qcMesh::Mesh MakeMolecularMesh(const Structure&, const qcMesh::MeshParams&);

} // namespace qchem
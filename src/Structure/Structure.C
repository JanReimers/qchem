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

    //! \brief Polymorphic deep copy that PRESERVES the concrete geometry (Atom / Molecule / UnitCell).  The
    //! production facade (qchem::MoleculeCalculation) clones the caller's structure through this, so a periodic
    //! UnitCell keeps its periodicity -- its Ewald ion-ion and uniform integration mesh -- instead of being
    //! sliced down to a finite Molecule.  Load-bearing for GPW (a real-space basis on a lattice).
    virtual std::shared_ptr<Structure> Clone() const=0;

    virtual std::string ID          () const;
    virtual size_t GetNumAtoms      () const=0;
    virtual int    GetNuclearCharge () const;
    virtual double GetNumElectrons  () const;
    virtual double GetNetCharge     () const;

    //! \brief The ion-ion (nuclear-nuclear) electrostatic energy, with the ION charge of each site given by \a zionOf
    //! (its true Z -> its core charge: identity all-electron, Z -> Zion for a pseudopotential).  Each geometry answers
    //! with its own sum: the direct Coulomb pair sum \f$\tfrac12\sum_{a\ne b}q_aq_b/|r_a-r_b|\f$ here (a finite site
    //! set: Atom, Molecule), the Ewald lattice sum on a periodic \c UnitCell (V-STRUCTURE-ISP F3: was a free
    //! \c NuclearRepulsion that asked \c isFinite() and down-cast to the concrete UnitCell).
    virtual double IonIonEnergy     (const std::function<double(int Z)>& zionOf) const;

    //! \brief The G = 0 NEUTRALISING-BACKGROUND term of a per-species form factor \a f: \f$\frac1\Omega\sum_a f(Z_a)\f$
    //! for a periodic structure (the per-volume mean a reciprocal-space term drops from its matrix and carries
    //! separately), exactly 0 for a finite one (no background, no dropped G = 0).  Asked by the PP G = 0 alignment
    //! (\c PWTerms_PP) and the GPW V_short cell mean (\c GPW_Evaluator); replaces the old \c isFinite() switch +
    //! \c SumFormFactors pair (V-STRUCTURE-ISP F3).
    virtual double G0Background     (const std::function<double(int Z)>& f) const {return 0.0;}

    //! \brief Does this structure's symmetry group contain the LATTICE TRANSLATIONS T?  The Vxc admission rule
    //! (\c Hamiltonian::AdmitVxc, doc/Pins.md pin 30) consumes it: a representation that imposes T (plane waves, a
    //! uniform cell mesh) is refused on a structure whose group lacks them.  False for a finite site set.
    //! \note the one structure-kind question F3 leaves on the site set; when the factory is handed the imposed
    //! group itself (V-STRUCTURE-ISP F2 / V-BS-SYMMETRY), this becomes that group's question.
    virtual bool   HasTranslationSymmetry() const {return false;}

    //! \brief Visit every atomic site as \c (Z,R,spinFlip) -- the assembly-level question seed and
    //! magnetic-decoration code ask (V1.19).  A consumer that needs "which species, where, which magnetic
    //! sublattice" states exactly that and stays off the concrete \c Atom's public fields; the storage
    //! stays here.  Sites are visited in atom-index order (consumers build per-atom-parallel tables).
    void ForEachSite(const std::function<void(int Z, const rvec3_t& R, bool spinFlip)>& f) const;

    //! \brief The real-space integration mesh natural to this geometry, at resolution \a mp.  The mesh TYPE
    //! is the STRUCTURE's responsibility, not the basis's: an Atom builds its single-centre (radial x angular)
    //! grid, a Molecule the multi-centre Becke grid, a periodic lattice a uniform / unit-cell-Becke grid.
    //! Pure virtual -- each geometry owns its most efficient mesh (no central dispatch).  A real-space
    //! orbital basis integrates a field on this via qcMesh::MatrixOverlap; plane waves own their G-grid and
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

    virtual std::shared_ptr<Structure> Clone() const override {return std::make_shared<Atom>(*this);}
    virtual size_t GetNumAtoms      () const override {return 1;}

    //! A single atom's integration mesh: the atom-centred radial x angular grid (no Becke partitioning).
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    virtual std::string   ID     () const override;
    virtual std::ostream& Write  (std::ostream&) const override;

    int     itsZ;      //Atomic number.
    double  itsCharge; //Net charge. Z-numElectrons.
    rvec3_t itsR;      //Spatial position.
    //! Collinear magnetic-configuration flip bit (doc/SCFSeedingPlan.md §10): the spin-SAD seed ASSEMBLY
    //! swaps this site's stored (majority, minority) channel pair -- the -m sublattice of an AFM pattern
    //! (e.g. MnO AFM-II).  Assembly-time data riding the structure (the MakeValenceStructure clone-and-
    //! mutate pattern); everything spin-agnostic ignores it.  Non-collinear (an SU(2) rotation per site)
    //! is the §4 generalization -- this bit is its collinear collapse.
    bool    itsSpinFlip = false;

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
    virtual std::shared_ptr<Structure> Clone() const override {return std::make_shared<Molecule>(*this);}
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
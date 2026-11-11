// File: Structure.C  A structure of atoms, Molecule, unit cell, solid all with spatial struture defined.
module;
#include <vector>
#include <memory>
#include <functional>
#include <cstddef>

export module qchem.Structure;
import qchem.Types;
import qchem.Streamable;
import qchem.Mesh;        // qcMesh::Mesh / MeshParams -- the real-space integration mesh a geometry produces

namespace qchem {

export class Structure;

//----------------------------------------------------------------------------------------------
//! \brief One site, a plain VALUE: the storage unit of every site set and what a client reads (V-STRUCTURE-ISP F1).
//! \a Z is the TRUE species (a pseudopotential maps it to Zion itself, never by editing Z); \a R the Cartesian
//! position (a.u.); \a netCharge = Z - (electrons the site carries), 0 for a neutral atom (GetNumElectrons /
//! GetNetCharge sum it); \a spinFlip the collinear magnetic-configuration flip bit: the spin-SAD seed ASSEMBLY swaps
//! this site's stored (majority, minority) channel pair -- the -m sublattice of an AFM pattern (e.g. MnO AFM-II,
//! doc/SCFSeedingPlan.md §10); everything spin-agnostic ignores it.  Aggregate, so `Site{Z, R}` builds one.
export struct Site
{
    int     Z;
    rvec3_t R         = {0,0,0};
    double  netCharge = 0.0;
    bool    spinFlip  = false;
};

//! \brief A read-only range over a structure's sites, for range-for: `for (const Site& s : st.Sites())`.  A VIEW: it
//! must not outlive the structure.  Sites come in atom-index order (consumers build per-atom-parallel tables).
export class SiteRange
{
public:
    class iterator   // a full forward iterator, so range-for AND the C++20 <ranges> adaptors compose on Sites()
    {
    public:
        using value_type      = Site;
        using difference_type = std::ptrdiff_t;
        iterator() = default;
        const Site& operator* () const;
        iterator&   operator++()                        {++itsI; return *this;}
        iterator    operator++(int)                     {iterator t=*this; ++itsI; return t;}
        bool        operator!=(const iterator& b) const {return itsI!=b.itsI;}
        bool        operator==(const iterator& b) const {return itsI==b.itsI;}
    private:
        friend class SiteRange;
        iterator(const Structure* s, size_t i) : itsS(s), itsI(i) {}
        const Structure* itsS=nullptr;
        size_t           itsI=0;
    };
    iterator    begin() const {return iterator(itsS,0);}
    iterator    end  () const {return iterator(itsS,size());}
    size_t      size () const;
    const Site& operator[](size_t i) const;           //!< the i-th site, atom-index order
private:
    friend class Structure;
    explicit SiteRange(const Structure& s) : itsS(&s) {}
    static const Site& At(const Structure* s, size_t i);   //!< the friend door to Structure::GetSite
    const Structure* itsS;
};

//----------------------------------------------------------------------------------------------
//! \brief THE SITE SET: an abstract interface over Atom (one site), Molecule and UnitCell (V-STRUCTURE-ISP F1).
//! Clients ASK: the sites (read-only, by value semantics -- no Atom* leaves this face), counts and charges, the
//! geometry's own integration mesh, its own ion-ion energy and G = 0 background (F3).  On a NON-const clone they may
//! MOVE the sites or SET their net charges (the clone-and-mutate pattern of the facades and setups).  Periodicity is
//! not here: the box is the UnitCell's, the translations the Lattice's (F2).
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
    //! \brief The FORCE of that energy, \f$-\partial E_{nn}/\partial R_a\f$ per site (GeometryOptimizationPlan §13.2), with the same charges
    //! and the same sum: the pair sum's closed form here, the Ewald sum's on a periodic \c UnitCell.
    virtual rvec3vec_t IonIonForces (const std::function<double(int Z)>& zionOf) const;
    //! \brief Its STRAIN derivative beside it (doc/Records/CellHalf_Design.md §5; user ruling §6.2, 2026-11-10, "the derivative version of an
    //! existing method, beside it"): \f$V_{ab}=\partial E_{nn}/\partial\varepsilon_{ab}\f$ under \f$R\to(1+\varepsilon)R\f$ (and a cell's
    //! \f$A\to(1+\varepsilon)A\f$), a GRADIENT in Ha, so a cell's stress is \f$V/\Omega\f$.  The pair sum's closed form here (well defined
    //! for a finite site set, but with no cell it is no stress); the Ewald sum's on a periodic \c UnitCell.
    virtual rmat3d_t   IonIonVirial (const std::function<double(int Z)>& zionOf) const;

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

    //! \brief The sites, read-only, in atom-index order (user ruling Q2: a by-value Site range).
    SiteRange Sites() const {return SiteRange(*this);}

    //! \brief Visit every atomic site as \c (Z,R,spinFlip) -- the callback form of \c Sites(), kept for its
    //! existing callers (V1.19; user ruling Q2).  Sites are visited in atom-index order.
    void ForEachSite(const std::function<void(int Z, const rvec3_t& R, bool spinFlip)>& f) const;

    //! \brief The real-space integration mesh natural to this geometry, at resolution \a mp.  The mesh TYPE
    //! is the STRUCTURE's responsibility, not the basis's: an Atom builds its single-centre (radial x angular)
    //! grid, a Molecule the multi-centre Becke grid, a periodic lattice a uniform / unit-cell-Becke grid.
    //! Pure virtual -- each geometry owns its most efficient mesh (no central dispatch).  A real-space
    //! orbital basis integrates a field on this via qcMesh::MatrixOverlap; plane waves own their G-grid and
    //! do not use it.
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const=0;

    virtual size_t GetAtomIndex(const rvec3_t&, double tol=0.0) const;

    //=== Mutators on a NON-const structure (clone-and-mutate); they act on the stored Sites ==========
    //! \brief Move every site to \a R (Cartesian a.u., one per site, atom-index order) -- a geometry step or a
    //! rigid shift on a cloned structure.
    void MoveSites(const std::vector<rvec3_t>& R);
    //! \brief Set each site's net charge from its species, \a netChargeOf(Z) -- a PP run's point structure carries
    //! net charge Z - Zion (MoleculeSetup::PointStructure).
    void SetNetCharges(const std::function<double(int Z)>& netChargeOf);

protected:
    friend class SiteRange;
    //! The storage primitives every concrete site set provides (a Molecule's vector, an Atom's one site).
    virtual const Site& GetSite(size_t i) const=0;
    virtual       Site& GetSite(size_t i)      =0;    //!< the mutators' door; never public
};

//! \brief The ONE-SITE site set: an isolated atom (AtomCalculation's structure).  It holds one Site; a Molecule
//! stores Sites, not Atoms (V-STRUCTURE-ISP F1: the old Atom was both a public-field record and a Structure).
export class Atom
    : public Structure
{
public:
    Atom(int Z);
    Atom(int Z, double netCharge);
    Atom(int Z, const rvec3_t& R);
    Atom(int Z, double netCharge, const rvec3_t& R);
    explicit Atom(const Site& s);
    //  Copy/move/assign are all defaulted: an Atom is a plain value type.

    virtual std::shared_ptr<Structure> Clone() const override {return std::make_shared<Atom>(*this);}
    virtual size_t GetNumAtoms      () const override {return 1;}

    //! A single atom's integration mesh: the atom-centred radial x angular grid (no Becke partitioning).
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    virtual std::string   ID     () const override;
    virtual std::ostream& Write  (std::ostream&) const override;

protected:
    virtual const Site& GetSite(size_t) const override {return itsSite;}
    virtual       Site& GetSite(size_t)       override {return itsSite;}
private:
    Site itsSite;
};

//! \brief A finite set of sites, stored BY VALUE (rule of zero: the vector copies, moves and destroys itself).
export class Molecule
    : public virtual Structure
{
public:
    Molecule() {};
    Molecule(const Structure& m);                     //!< copies \a m's sites (any site set)
    virtual std::shared_ptr<Structure> Clone() const override {return std::make_shared<Molecule>(*this);}
    //! \brief Add a site (was Insert(Atom*) taking ownership of a raw new).
    virtual void   Insert     (const Site& s);
    virtual size_t GetNumAtoms() const override {return itsSites.size();}

    //! A molecule's integration mesh: the multi-centre Becke-partitioned grid.
    virtual qcMesh::Mesh CreateIntegrationMesh(const qcMesh::MeshParams&) const override;

    virtual std::ostream& Write(std::ostream&) const override;

protected:
    virtual const Site& GetSite(size_t i) const override {return itsSites[i];}
    virtual       Site& GetSite(size_t i)       override {return itsSites[i];}

private:
    std::vector<Site> itsSites;
};

inline const Site& SiteRange::At(const Structure* s, size_t i) {return s->GetSite(i);}
inline const Site& SiteRange::iterator::operator*() const      {return At(itsS,itsI);}
inline size_t      SiteRange::size() const                     {return itsS->GetNumAtoms();}
inline const Site& SiteRange::operator[](size_t i) const       {return At(itsS,i);}

// One site's ID / console line -- the old Atom::ID / Atom::Write, shared by Atom and the multi-site Write/ID.
// Module-INTERNAL (not exported).
std::string   SiteID   (const Site&);
std::ostream& WriteSite(std::ostream&, const Site&);

// Becke fuzzy-Voronoi integration mesh for a finite structure -- the shared implementation behind
// Atom::/Molecule::CreateIntegrationMesh.  Module-INTERNAL (NOT exported): callers use CreateIntegrationMesh,
// the public per-geometry entry point; this is the algorithm both geometries delegate to.  See Imp/Molecule.C
// (Becke 1988; coincident atoms R_ab=0 -> mu=0 so a coincident dimer integrates to the single-atom result).
qcMesh::Mesh MakeMolecularMesh(const Structure&, const qcMesh::MeshParams&);

} // namespace qchem
// File: Structure.C  A structure of atoms, Molecule, unit cell, solid all with spatial struture defined.
module;
#include <vector>
#include <memory>
#include <functional>
#include <cstddef>

export module qchem.Structure;
import qchem.Types;
import qchem.Streamable;
import qchem.Symmetry.Group;   // Symmetry::GroupPtr / GroupTolerance (MakeSymmetry); a client that asks the group imports it itself

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
        iterator&   operator++();
        iterator    operator++(int);
        bool        operator!=(const iterator& b) const;
        bool        operator==(const iterator& b) const;
    private:
        friend class SiteRange;
        iterator(const Structure* s, size_t i);
        const Structure* itsS=nullptr;
        size_t           itsI=0;
    };
    iterator    begin() const;
    iterator    end  () const;
    size_t      size () const;
    const Site& operator[](size_t i) const;           //!< the i-th site, atom-index order
private:
    friend class Structure;
    explicit SiteRange(const Structure& s);
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
    virtual ~Structure();

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
    virtual double G0Background     (const std::function<double(int Z)>& f) const;

    //! \brief The structure's own symmetry group G_struct, in the Cartesian frame (lane H Group face, migration step 2), detected at
    //! \a tol.  MAKE, not Get: it computes, and the result depends on \a tol (a later cache must key on it).  Each kind issues its own:
    //! an \c Atom O(3) about its site (the atom path imposes it by spherically averaging its open shell); a \c Molecule the ABELIAN
    //! subgroup of its point group, the largest we can enumerate today and the group the molecular path imposes (user, 2026-11-19;
    //! enumerating the full point group is an OpenWork [G3] row); a \c UnitCell its space group, lattice and every coset included.
    //! The group is SPATIAL (Space::R): a site's spinFlip does not enter (the magnetic factor is migration step 8).
    virtual Symmetry::GroupPtr MakeSymmetry(const Symmetry::GroupTolerance& tol) const=0;

    //! \brief The sites, read-only, in atom-index order (user ruling Q2: a by-value Site range).
    SiteRange Sites() const;

    //! \brief Visit every atomic site as \c (Z,R,spinFlip) -- the callback form of \c Sites(), kept for its
    //! existing callers (V1.19; user ruling Q2).  Sites are visited in atom-index order.
    void ForEachSite(const std::function<void(int Z, const rvec3_t& R, bool spinFlip)>& f) const;

    virtual size_t GetAtomIndex(const rvec3_t&, double tol=0.0) const;

    //=== Mutators on a NON-const structure (clone-and-mutate); they act on the stored Sites ==========
    //! \brief Move every site to \a R (Cartesian a.u., one per site, atom-index order) -- a geometry step or a
    //! rigid shift on a cloned structure.
    void MoveSites(const std::vector<rvec3_t>& R);
    //! \brief ADD \a extraChargeOf(Z) to each site's net charge, keeping the charge the site already states -- a PP run's point
    //! structure adds the core charge Z - Zion (MoleculeSetup::PointStructure), so a stated ion (OH-) stays an ion.  It once REPLACED
    //! the stated charge, and OH- with GTH ran neutral (lane G, 2026-11-21).
    void AddNetCharges(const std::function<double(int Z)>& extraChargeOf);

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

    virtual std::shared_ptr<Structure> Clone() const override;
    virtual size_t GetNumAtoms      () const override;
    //! O(3) about the site: the atom path imposes it (its open shell is spherically averaged, so its density is spherical).
    virtual Symmetry::GroupPtr MakeSymmetry(const Symmetry::GroupTolerance& tol) const override;

    virtual std::string   ID     () const override;
    virtual std::ostream& Write  (std::ostream&) const override;

protected:
    virtual const Site& GetSite(size_t) const override;
    virtual       Site& GetSite(size_t)       override;
private:
    Site itsSite;
};

//! \brief A finite set of sites, stored BY VALUE (rule of zero: the vector copies, moves and destroys itself).
export class Molecule
    : public virtual Structure
{
public:
    Molecule();
    //! Copies \a m's sites (any site set).  EXPLICIT (USER RULING 2026-11-23): implicit, it let an abstract \c Structure -- a cell, an
    //! atom -- convert silently into a molecule wherever a \c const \c Molecule& was taken (16 mesh calls did, V-MESH-FACTORY step 4b).
    explicit Molecule(const Structure& m);
    virtual std::shared_ptr<Structure> Clone() const override;
    //! \brief Add a site (was Insert(Atom*) taking ownership of a raw new).
    virtual void   Insert     (const Site& s);
    virtual size_t GetNumAtoms() const override;
    //! The abelian subgroup of the point group, about the centroid (see the base's declaration).  THROWS for an empty site set,
    //! whose group would be every isometry: there is no finite group to issue.
    virtual Symmetry::GroupPtr MakeSymmetry(const Symmetry::GroupTolerance& tol) const override;

    virtual std::ostream& Write(std::ostream&) const override;

protected:
    virtual const Site& GetSite(size_t i) const override;
    virtual       Site& GetSite(size_t i)       override;

private:
    std::vector<Site> itsSites;
};

// One site's ID / console line -- the old Atom::ID / Atom::Write, shared by Atom and the multi-site Write/ID.
// Module-INTERNAL (not exported).
std::string   SiteID   (const Site&);
std::ostream& WriteSite(std::ostream&, const Site&);

} // namespace qchem
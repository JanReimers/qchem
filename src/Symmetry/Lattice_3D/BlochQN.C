// File: Symmetry/Lattice_3D/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iosfwd>
#include <string>
#include <vector>
export module qchem.Symmetry.Lattice_3D.BlochQN;
export import qchem.Types;
export import qchem.Symmetry;
export import qchem.Symmetry.SelectionRule;   // MeshShift IS the lattice's selection rule (k -> k+q)
export import qchem.Outcome;   // CommensurateShifts: an incommensurate q-mesh is a configuration error, not a throw
//---------------------------------------------------------------------------------
//
//  Translational symmetry, Bloch function wave vector.
//

namespace qchem::Symmetry::Lattice_3D {

export class BlochQN;

//! \brief A wave-vector SHIFT q ON A k-MESH'S OWN GRID: an integer number of grid steps \f$\Delta ik\f$ on a
//! mesh of \f$N\f$ divisions, \f$q=\Delta ik/N\f$ (fractional reciprocal coordinates).
//! (doc/LinearResponsePlan.md §3 row S1, stage R0.)
//!
//! ★ WHY IT IS A TYPE AND NOT A \c rvec3_t.  A linear response couples block k to block k+q, and the k+q
//! block must EXIST -- which it does exactly when q is a difference of mesh points.  A float q can be off the
//! mesh; this type cannot, because only a mesh point (\c BlochQN::CommensurateShifts) can construct one.  An
//! off-mesh q is therefore UNREPRESENTABLE rather than checked (build failure over runtime failure, user
//! 2026-08).  For every MeshShift, k+q IS a mesh point, shifted Monkhorst-Pack meshes included, because q is
//! a difference of two of them.
//! \note In our Bloch gauge (phase \f$e^{i\mathbf k\cdot\mathbf R_n}\f$ over integer lattice offsets,
//! GPW_Evaluator / LatticeSum1E) \f$\phi_{k+G}\equiv\phi_k\f$, so "k+q modulo a reciprocal lattice vector"
//! is a pure index map with no G-phase to carry -- unlike QE's \f$e^{i(\mathbf k+\mathbf G)\cdot\mathbf r}\f$
//! basis, which needs its \c ikqs table.
export class MeshShift : public virtual qchem::Symmetry::SelectionRule
{
public:
    //! The lattice selection rule: \a bra is \f$k+q\f$ for \a ket (both Bloch points of THIS mesh).
    virtual bool Couples(const qchem::Symmetry::Symmetry& bra, const qchem::Symmetry::Symmetry& ket) const override;
    ivec3_t Grid () const {return N;}    //!< the k-mesh divisions this shift lives on
    ivec3_t Steps() const {return d;}    //!< \f$\Delta ik\f$, reduced into [0, N)
    rvec3_t q    () const;               //!< \f$\Delta ik/N\f$, fractional reciprocal coordinates
    bool    IsZero() const {return d.x==0 && d.y==0 && d.z==0;}
private:
    friend class BlochQN;
    MeshShift(ivec3_t _N, ivec3_t _d);
    ivec3_t N, d;
};
export std::ostream& operator<<(std::ostream&, const MeshShift&);

export class BlochQN : public virtual qchem::Symmetry::Symmetry
{
public:
    //! \a _N = BZ-grid divisions, \a _ik = integer grid index (0..N, used for the sequence index / identity),
    //! \a _shift = the fractional Monkhorst-Pack offset in grid steps so \f$k=(ik+shift)/N\f$: \a shift=0 is the
    //! Γ-centred grid; \a shift=½ is the classic MP offset (\f$k=\pm¼\f$ for \f$N=2\f$, i.e. CP2K's default).
    //! \a _weight = the block's BZ integration weight \f$w_k\f$ (the k-mesh layer's currency: \f$1/N\f$ on an
    //! unfolded mesh, star/\f$N\f$ for an IBZ wedge representative; \f$\sum_\mathrm{blocks} w_k=1\f$).
    //! INTERNALLY the ATOM SHELL CONVENTION applies (user design, 2026-08-04): the star multiplicity
    //! \f$w_k N_\mathrm{mesh}\f$ (asserted integer) is carried as the block's spatial DEGENERACY -- one stored
    //! representative stands for its whole k-star, exactly like an atom's l-shell stores one radial with
    //! degeneracy \f$2l+1\f$ -- and \c GetWeight returns the UNIFORM per-point \f$1/N_\mathrm{mesh}\f$.  The
    //! physical invariant \f$w_k\times\f$(per-block quantity) is unchanged: occupations/densities/entropy pick
    //! up the star through the degeneracy-driven fill (\c Crystal_EC::GetN scales with it), the weight drops to
    //! \f$1/N\f$, and an unfolded mesh (star = 1) is bit-identical to the old convention.
    BlochQN(ivec3_t _N, ivec3_t _ik, double _weight=1.0, rvec3_t _shift={0,0,0});
    virtual size_t SequenceIndex() const;
    //! Spatial degeneracy = the k-star multiplicity of this block (1 on an unfolded mesh; the star size on
    //! an IBZ wedge).  Spin rides on top via \c Irrep, so a wedge band level holds \f$2\cdot\f$star electrons.
    virtual size_t GetDegeneracy() const {return star;}
    virtual size_t GetPrincipleOffset() const  {return 1;}
    virtual double GetWeight() const; // the UNIFORM per-point sampling weight 1/N_mesh (star lives in GetDegeneracy)
    //! A k-STAR of band levels is one symmetry-degenerate shell of the crystal group -- let the
    //! EnergyLevels reporting layer merge equal-eigenvalue levels across k-blocks (Symmetry doc).
    virtual bool   MergeAcrossIrreps() const {return true;}
    //! \brief TRIM test: real iff \f$2k\equiv 0\pmod{\mathrm{recip.\ lattice}}\f$, i.e.
    //! \f$N_i\,|\,2(ik_i+\mathrm{shift}_i)\f$ per component -- EXACT integer arithmetic on the ctor's
    //! integer grid data, no float-k tolerance (doc/RealComplexPlan.md Step 1; the phases themselves are
    //! exactly \f$\pm1\f$ there by Step 0).  \f$\Gamma\f$ and zone-boundary half-integer k answer yes; a
    //! shifted (MP \f$\mathrm{shift}=\tfrac12\f$, even N) mesh point does not.
    virtual bool   IsReal() const {return isReal;}
    virtual std::ostream&  Write(std::ostream&) const;

    rvec3_t   Getk() const {return k;}
    //! \brief The q-mesh of \a Nq divisions ON THIS k-mesh: every \f$q=iq/N_q\f$, as MeshShifts.  Its
    //! \f$\prod N_q\f$ points are the supercell images a monochromatic response sums over.
    //! FAILS (a value, not a throw -- the caller chose the q-mesh and can choose another) unless \f$N_q\f$
    //! divides \f$N\f$ per axis: that is the commensurability every DFPT code requires (QE's hp.x: "limited to
    //! q point grids that are commensurate with the k point grid").
    Outcome<std::vector<MeshShift>,std::string> CommensurateShifts(ivec3_t Nq) const;
    //! Is THIS point \f$k+q\f$ for \a k (same mesh, same Monkhorst-Pack shift), modulo a reciprocal lattice
    //! vector?  Exact integer arithmetic; a \a q or \a k from a different mesh answers false.
    bool IsShiftOf(const BlochQN& k, const MeshShift& q) const;

private:
    ivec3_t N;      //This is the Brillouin zone grid size which gives context for the k vector. Used for calculating the sequence index.
    ivec3_t ik;     //Integer rep. of k.
    rvec3_t k;      //Real values.
    rvec3_t shift;  //The Monkhorst-Pack offset in grid steps (0 = Γ-centred) -- two points share a mesh only if equal.
    size_t  star;   //k-star multiplicity w_k·N_mesh (the block's spatial degeneracy; 1 on an unfolded mesh).
    bool    isReal; //TRIM fact N_i | 2(ik_i+shift_i), computed EXACTLY in the ctor (see IsReal).
};

//! \brief Pry the Bloch wave vector \f$k\f$ (fractional) out of an abstract symmetry handle.
//! Throws std::bad_cast if the symmetry is not a BlochQN.  Mirrors Symmetry::Atom::Getl / Getκ:
//! an IBS constructor is handed an abstract \c sym_t and uses this helper to extract the one
//! piece of concrete information it needs (here, the crystal momentum).
export rvec3_t Getk(const sym_t&);
export rvec3_t Getk(const qchem::Symmetry::Symmetry&);
//! Pry the same way (throws std::bad_cast on a non-Bloch handle): is \a kq the point \f$k+q\f$ for \a k?
export bool IsShiftOf(const qchem::Symmetry::Symmetry& kq, const qchem::Symmetry::Symmetry& k, const MeshShift& q);
//! The q-mesh of \a Nq divisions on \a anyK's k-mesh (see \c BlochQN::CommensurateShifts).
export Outcome<std::vector<MeshShift>,std::string> CommensurateShifts(const qchem::Symmetry::Symmetry& anyK, ivec3_t Nq);
} // namespace qchem::Symmetry::Lattice_3D


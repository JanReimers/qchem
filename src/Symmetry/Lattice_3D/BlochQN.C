// File: Symmetry/Lattice_3D/BlochQN.C  A Quantum Number translational symmetry, i.e. a wave vector.
module;
#include <iosfwd>
export module qchem.Symmetry.Lattice_3D.BlochQN;
export import qchem.Types;
export import qchem.Symmetry;
//---------------------------------------------------------------------------------
//
//  Translational symmetry, Bloch function wave vector.
//

namespace qchem::Symmetry::Lattice_3D {
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
    virtual std::ostream&  Write(std::ostream&) const;

    rvec3_t   Getk() const {return k;}

private:
    ivec3_t N;      //This is the Brillouin zone grid size which gives context for the k vector. Used for calculating the sequence index.
    ivec3_t ik;     //Integer rep. of k.
    rvec3_t k;      //Real values.
    size_t  star;   //k-star multiplicity w_k·N_mesh (the block's spatial degeneracy; 1 on an unfolded mesh).
};

//! \brief Pry the Bloch wave vector \f$k\f$ (fractional) out of an abstract symmetry handle.
//! Throws std::bad_cast if the symmetry is not a BlochQN.  Mirrors Symmetry::Atom::Getl / Getκ:
//! an IBS constructor is handed an abstract \c sym_t and uses this helper to extract the one
//! piece of concrete information it needs (here, the crystal momentum).
export rvec3_t Getk(const sym_t&);
export rvec3_t Getk(const qchem::Symmetry::Symmetry&);
} // namespace qchem::Symmetry::Lattice_3D


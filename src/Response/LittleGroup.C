// File: Response/LittleGroup.C  The k-mesh of ONE perturbation reduced by its little group (doc/KAxisPlan.md P5, D5-a).
//
// A response to the perturbation V_J e^{iq.R} needs dD on every (k+q, k) pair of the full mesh -- but an op g of the ground
// state's group that maps q to itself (mod G) and the perturbation onto itself (up to a phase chi_g) maps the solution onto
// itself:
//
//     dD_{gk} = chi_g U_g(k+q) dD_k U_g(k)^H,      U_g(k) = C_{gk}^H S_{gk} M_g(k) C_k   (the op on the ORBITAL basis)
//
// with M_g(k) the basis's OrbitalImage (P4).  So the unknowns are one block per ORBIT of that little group, and the rest is
// unfolded.  Two values here:
//   * OrbitalImages -- the ground state's group as MEASURED on the reference (KAxisPlan §5 ruling 1): an op is kept only if
//     U_g maps orbitals onto orbitals of the same eigenvalue AND occupation (max|U(eps'_i - eps_j)|, max|U(f'_i - f_j)|) and
//     is unitary, each to the tolerance.  That is exactly what the response needs: R0 then commutes with g.  One per run.
//   * LittleGroupReduction -- for one (q, perturbation): the ops that also keep q and V (the phase chi_g MEASURED, never
//     derived), the orbits, and the reduced PACKING the Krylov solver works on.  Each representative carries sqrt(orbit
//     size), so the packing is an ISOMETRY of the symmetric subspace of the full one: GMRES takes the same steps as on the
//     full mesh (in exact arithmetic), and its residual means the same thing.  A representative with a non-trivial
//     stabilizer is PROJECTED onto the stabilizer-invariant subspace on every pack (rounding would otherwise drift out).
//
// ⚠ COLLINEAR ONLY: a spin-flip op (sigma) is applied as a relabel of the up/down BLOCKS, with the spatial U_g.  Non-collinear
// magnetism needs the SU(2) rep on a spinor block (M (x) U(g)), which changes the image construction (MakeOrbitalImages) and
// this unfolding.
module;
#include <iosfwd>
#include <vector>
export module qchem.Response.LittleGroup;
export import qchem.Response.Reference;
export import qchem.Symmetry.SymOp;

export namespace qchem::Response
{

//! One op on the reference's blocks, in the orbital basis: block \c b goes to block \c target[b], its orbitals by
//! \c U[b] (target orbitals x source orbitals).
struct OrbitalOp
{
    Symmetry::SymOp     op;
    std::vector<size_t> target;
    std::vector<cmat_t> U;
};

//! The ground state's group as measured on the reference (ruling 1): the ops kept, and what each rejected op missed by.
struct OrbitalImages
{
    std::vector<OrbitalOp> ops;           //!< the KEPT ops (the identity among them)
    size_t                 tried=0;       //!< ops offered
    double                 worstKept=0;   //!< the largest residual among the kept ops
    //! One rejected op: what it missed by, on which test and source block (inf = an image block that is not stored, or an op the
    //! basis cannot image) -- named so a systematic split (e.g. every op with a fractional translation) is readable off the log.
    struct Rejection { Symmetry::SymOp op; double resid; const char* test; size_t block; };
    std::vector<Rejection> rejected;
    double                 tol=1e-6;
    std::ostream& Write(std::ostream&) const;
};

//! \brief The little group of ONE perturbation \a V (on \a rule's pairs, orbital basis) and the k-blocks it reduces to.
//! THROWS (std::logic_error) if the kept ops are not a group -- the measured group of a symmetric state always is.
class LittleGroupReduction
{
public:
    LittleGroupReduction(const Reference& ref, const OrbitalImages& images, const SelectionRule& rule, const BlockPairs& V,
                         double tol=1e-6);
    size_t NumOps() const {return itsOps.size();}
    //! The orbit representatives (ket block indices, in reference order) -- the blocks a reduced kernel must answer.
    const std::vector<size_t>& Representatives() const {return itsReps;}
    //! The reduced vector of a symmetric full \a x: each representative projected, scaled by sqrt(orbit), column-major.
    cvec_t     Pack  (const BlockPairs& x) const;
    //! The inverse of Pack, unfolded onto every block.
    BlockPairs Unpack(const cvec_t& v) const;
    //! Every block from the representatives' entries of \a onReps (the others are ignored): \f$x_{gk}=\chi_gU x_kU^\dagger\f$.
    BlockPairs Unfold(const BlockPairs& onReps) const;
    std::ostream& Write(std::ostream&) const;

private:
    cmat_t Apply(size_t op, size_t ket, const cmat_t& x) const;   //!< chi U(bra) x U(ket)^H for the op's image of the pair
    cmat_t Project(size_t rep, const cmat_t& x) const;            //!< the stabilizer average of x

    const Reference&          itsRef;
    std::vector<size_t>       itsPartner;   //!< per ket block: its bra (the rule's)
    std::vector<OrbitalOp>    itsOps;       //!< the little group
    std::vector<dcmplx>       itsChi;       //!< per op: the perturbation's measured phase
    std::vector<size_t>       itsReps;
    std::vector<size_t>       itsRepOf, itsOpOf;          //!< per block: its representative and the op reaching it
    std::vector<std::vector<size_t>> itsStab;             //!< per representative (indexed like itsReps): its stabilizer
    std::vector<double>       itsScale;                   //!< per representative: sqrt(orbit size)
    size_t                    itsDim=0;
    double                    itsWorstV=0;                //!< the worst perturbation residual among the kept ops
};

} // namespace

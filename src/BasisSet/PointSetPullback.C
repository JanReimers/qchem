// File: BasisSet/PointSetPullback.C  A block's point-set quadrature pullback: the capability an atom-centred XC quadrature asks of an orbital
// block to differentiate its matrix (lane B 2026-11-19, READY 4: the δ XC force on the Becke mesh).
module;
#include <vector>
export module qchem.BasisSet.PointSetPullback;
import qchem.Structure;                       // the sites the functions sit on
export import qchem.Types;                    // rvec3vec_t, rmat3d_t, hmat_t<T>
export import qchem.AD;                       // NoGradient: a block with no such derivative refuses
export import qchem.Outcome;
export import qchem.BasisSet.KFoldable;       // KFoldTerm (the k-fold of the pullback)

export namespace qchem::BasisSet
{

//! \brief THE DERIVATIVE OF A POINT-SET QUADRATURE over this block's functions: for \f$E=\sum_q c_q\,\mathrm{Re}\sum_{ij}D_{ji}
//! \overline{\phi_i(r_q)}\phi_j(r_q)\f$ (Re Tr(D M), M the matrix a quadrature with weights x field \a c assembles at \a pts) at FIXED \a c
//! and Hermitian \a D:  \c gradient per site of \a st with the functions' centres moving and the points held; \c moment
//! \f$\sum_C(\partial E/\partial C)\otimes(C-r_q)\f$ over every unfolded centre (origin-free); \c perPoint \f$\partial E/\partial r_q\f$ with the
//! centres held.  A quadrature whose points ride their owners combines them (\c LatticeCollocation::BlochPointPullback's doc).
//! CLIENT: the atom-centred XC sampler's MakeMatrixPullback / MakeMatrixVirial, by an abstract -> abstract cross-cast from the block it is
//! handed per call (as its \c Matrix cross-casts to the DFT tier).  A block without the face is a named MissingGradients entry.
template <class T> class PointSetPullback
{
public:
    virtual ~PointSetPullback() = default;
    struct Result { rvec3vec_t gradient; rmat3d_t moment; rvec3vec_t perPoint; };
    virtual Outcome<Result,AD::NoGradient> MakePointSetPullback(const rvec3vec_t& pts, const rvec_t& c, const hmat_t<T>& D, const Structure* st) const=0;
    //! Can \c MakePointSetPullback answer?  A block that implements it says yes; a DECORATOR that forwards it asks what it wraps (the SALC
    //! block over a raw AO basis).  CLIENT: the XC sampler's \c CanPullbackMatrix (the term's MissingGradients entry, at setup).
    virtual bool CanPointSetPullback() const {return true;}
    //! \brief THE k-FOLD (user ruling 2026-11-22 "yes"): the SUM of \c MakePointSetPullback over \a terms, blocks of MY fold domain (me
    //! among them; \c KFoldable), each with its own D -- the functions evaluated once at each point for all of them.  CLIENT: the XC
    //! sampler's multi-block pullback, which groups the blocks.  Default: a named refusal (a block that does not fold).
    virtual Outcome<Result,AD::NoGradient> MakePointSetPullbackFolded(const rvec3vec_t& /*pts*/, const rvec_t& /*c*/,
                                                                     const std::vector<KFoldTerm>& /*terms*/, const Structure* /*st*/) const
    {return Outcome<Result,AD::NoGradient>::Fail({"this block's point-set pullback does not fold over k"});}
};

} // namespace qchem::BasisSet

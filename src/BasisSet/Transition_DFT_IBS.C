// File: BasisSet/Transition_DFT_IBS.C  The DFT integrals of a (k+q, k) TRANSITION density -- the basis side of a
// q != 0 linear response (doc/LinearResponsePlan.md §3d, row B2).
//
// A transition density of wave vector q is
//     δρ(r) = Σ_ij δD_ij χ_i^{k+q}(r) conj(χ_j^k(r)),
// a Bloch-q function: NOT lattice-periodic, so neither the ground-state 3-centre tensors (Overlap3C /
// Repulsion3C, which assume a Hermitian D on ONE block) nor their G-space maps can carry it.  These faces are
// the pair-shaped siblings, and they are CROSS-CAST CAPABILITIES (the RealBlock idiom): the response terms ask a
// block for them and a basis without them fails loudly.  At q = 0 with a Hermitian δD on one block each reduces
// to its ground-state counterpart (gated in UTResponse).
//
// TWO FACES, because the two XC quadratures hold their tables in different places:
//  * \c Transition_DFT_IBS -- on the KET orbital block: the collocated routes (Hartree on the fit ball, and the
//    raw raster feed of the uniform-mesh XC).  The block owns the collocation machinery.
//  * \c Transition_Overlap3C -- on a POINT (δ) fit basis: the Becke/δ route, whose Φ tables that basis owns.
//
// NOTHING HERE IS FOLDED OR STAR-AVERAGED, whatever symmetry the run imposed: a perturbation breaks the group
// (§3d finding 5).  That is by construction -- these are new code with no fold in them.
module;
#include <string>
export module qchem.BasisSet.Transition_DFT_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;   // Orbital_DFT_IBS, cFIT_CD_ABS / cFIT_SF_ABS, ΔGq_Map
import qchem.Types;

export namespace qchem::BasisSet
{

//! \brief The (k+q, k) transition-density integrals of a COLLOCATING orbital block -- asked of the KET block.
//!
//! \a bra is the (k+q) block, \a dD is bra rows x ket columns and NOT Hermitian, and \a q is FRACTIONAL (the
//! selection rule's \c WaveVectorShift::q -- one representative for every pair).  The transition density is
//! collocated through the ket's AO basis, so bra and ket must be built over the SAME AO set (§3d finding 6):
//! every method THROWS when \c TransitionBasisID differs.  Periodic by construction, hence \c dcmplx throughout.
class Transition_DFT_IBS
{
public:
    virtual ~Transition_DFT_IBS() = default;
    //! The identity of the AO set my Bloch functions are built over.  A transition pair can be collocated
    //! through ONE basis only when bra and ket answer the same string (a per-k ortho drop would break it
    //! silently -- which is what the pin-22 vet trim prevents, and this checks).
    virtual std::string TransitionBasisID() const=0;
    //! \f$\delta V_H(G+q)=4\pi\,\delta\tilde\rho(G+q)/|G+q|^2\f$ on \a fit's ball, G = 0 KEPT iff q != 0 --
    //! the transition sibling of \c Repulsion3C(fit).apply.
    virtual ΔGq_Map TransitionRepulsion(const cFIT_CD_ABS& fit, const Transition_DFT_IBS& bra,
                                        const mat_t<dcmplx>& dD, const rvec3_t& q) const=0;
    //! Its EXACT adjoint: \f$h_{ij}=\langle\chi_i^{k+q}|\sum_G V(G+q)e^{i(G+q)\cdot r}|\chi_j^k\rangle\f$, bra x ket.
    virtual mat_t<dcmplx> TransitionPotential(const cFIT_CD_ABS& fit, const Transition_DFT_IBS& bra,
                                              const ΔGq_Map& V) const=0;
    //! The RAW feed on \a fit's raster: the PERIODIC PART \f$u=e^{-iq\cdot r}\delta\rho\f$ at the raster points --
    //! the transition sibling of \c Overlap3C(fit).applyRaw (collocated per ladder level, coarse levels
    //! spectrally transferred in).  Complex; real only at q = 0 with a Hermitian δD.
    virtual cvec_t TransitionOnGrid(const cFIT_SF_ABS& fit, const Transition_DFT_IBS& bra,
                                    const mat_t<dcmplx>& dD, const rvec3_t& q) const=0;
    //! Its EXACT adjoint for a periodic field \a v at the raster points:
    //! \f$h_{ij}=\langle\chi_i^{k+q}|e^{iq\cdot r}v(r)|\chi_j^k\rangle\f$, bra x ket.
    virtual mat_t<dcmplx> TransitionGridAdjoint(const cFIT_SF_ABS& fit, const Transition_DFT_IBS& bra,
                                                const cvec_t& v, const rvec3_t& q) const=0;
};

//! \brief The same pair for a POINT (δ) fit basis, which owns the orbital value tables \f$\Phi_{ai}=\chi_i(r_a)\f$
//! -- the Becke / δ-quadrature XC route.  No phase is needed here (§3d finding 4): at a mesh point
//! \f$\delta\rho(r_a)=[\Phi^{\rm bra}\delta D\,\Phi^{{\rm ket}\dagger}]_{aa}\f$ is the Bloch-q function itself, and the
//! adjoint integrates a lattice-periodic integrand over the cell's partition.
class Transition_Overlap3C
{
public:
    virtual ~Transition_Overlap3C() = default;
    //! \f$\delta\rho(r_a)=\sum_{ij}\Phi^{\rm bra}_{ai}\,\delta D_{ij}\,\overline{\Phi^{\rm ket}_{aj}}\f$ at my points.
    virtual cvec_t TransitionForward(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                     const mat_t<dcmplx>& dD) const=0;
    //! Its exact adjoint: \f$h_{ij}=\sum_a w_a\,\overline{\Phi^{\rm bra}_{ai}}\,v_a\,\Phi^{\rm ket}_{aj}\f$, bra x ket.
    virtual mat_t<dcmplx> TransitionAdjoint(const Orbital_DFT_IBS<dcmplx,dcmplx>& bra, const Orbital_DFT_IBS<dcmplx,dcmplx>& ket,
                                            const cvec_t& v) const=0;
};

} // namespace

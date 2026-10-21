// File: Hamiltonian/HubbardChannels.C  The +U term's PROJECTORS, as a public capability: what a linear
// response PERTURBS and MEASURES (doc/LinearResponsePlan.md §2, stage R0).
//
// WHY A FACE HERE.  A response probe perturbs Hubbard channel J with \f$\alpha\hat P_J\f$ and measures the
// channel occupations \f$n_I={\rm Tr}\,\hat P_I\rho\f$.  Both halves need the SAME projector the +U energy
// uses -- the Löwdin / atomic / ortho-atomic flavour the run chose -- or the response and the ground state
// disagree about what "the Ni 3d occupation" is (insight 4 of doc/HubbardUPlan.md: "if A7 ever disagrees
// with hp.x, the FIRST suspect must be a different projector").  The projector lives behind the .Internal.
// Hubbard module, and the response library sits ABOVE qcHamiltonian, so the Hamiltonian hands out this
// capability (tHamiltonian::GetHubbardChannels) and the term implements it: projector consistency by
// CONSTRUCTION, not by discipline.  The same shape as HubbardUEstimator: the data crosses the library
// boundary, the mechanism does not.
module;
#include <cstddef>
#include <string>
#include <vector>
export module qchem.Hamiltonian.HubbardChannels;
import qchem.BasisSet.Orbital_DFT_IBS;   // the block an orbital set lives on
import qchem.Types;

export namespace qchem::Hamiltonian
{

//! One Hubbard channel: a manifold, named the way a report quotes it.
struct HubbardChannel
{
    size_t site=0;   //!< the manifold's cell site (Structure order)
    int    l=2;      //!< its shell
    //! Does the manifold carry a +U potential (U != 0 in some slot)?  The rest are PROJECTOR SPECTATORS -- in
    //! the set only so the ortho-atomic projectors are orthogonalised against them (QE's ortho-atomic block).
    //! A linear-response U perturbs the carrying set by default: hp.x's "Hubbard sites" (LinearResponsePlan Q10).
    bool   carriesU=false;
};

//! One U SLOT of a manifold as a subspace of its projector functions: what an orbital-RESOLVED response perturbs and
//! measures (OpenWork §2 "Response: REPORT orbital-resolved channels", 2026-10-21).  The slots are the +U term's own --
//! the (site irrep, grey parent) isotypic components its \c Uirrep applies to -- so a resolved U maps onto a slot.
struct HubbardSlot
{
    std::string label;      //!< "Γk" (site irrep k), with the grey parent when the run is magnetic: "Γk<p"
    size_t      irrepDim=1; //!< the site irrep's dimension (the slot's dim is a multiple: the irrep's multiplicity)
    rmat_t      P;          //!< \f$m_M\times m_M\f$ orthogonal projector onto the slot, in the amplitudes' row basis
};
//! A manifold's slots, and the names of its \f$m_M\f$ projector functions (amplitude rows) so a report can say what
//! a slot IS ("z2 x2-y2" = e_g).  BOTH EMPTY when the term has no symmetry for the manifold (the population form):
//! the manifold is then one subspace.
struct HubbardSlots
{
    std::vector<std::string> components;
    std::vector<HubbardSlot> slots;
};

//! \brief The projector set of the run's +U term, one CHANNEL per Hubbard manifold, in the term's manifold
//! order.  \f$\hat P_J=\sum_\mu|w_\mu\rangle\langle w_\mu|\f$ over the manifold's projector functions
//! \f$w_\mu\f$; a channel is described to a client entirely by the orbitals' AMPLITUDES on those functions.
class HubbardChannels
{
public:
    virtual ~HubbardChannels() = default;
    virtual std::vector<HubbardChannel> Channels() const = 0;
    //! \brief Per channel, the amplitudes \f$\ell=\langle w_\mu|\psi_i\rangle=T_M^\dagger C\f$
    //! (\f$m_M\times n_{\rm orb}\f$) of the orbitals whose AO coefficients are the columns of \a C, on
    //! \a block.  Exactly what the +U occupation matrix is made of: \f$n_M=\sum f\,\ell\ell^\dagger\f$.
    //! \note For a Bloch block the functions are the Bloch sums over lattice translations (phase
    //! \f$e^{i\mathbf k\cdot\mathbf R_n}\f$), so \f$\ell\f$ is a per-cell amplitude.
    virtual std::vector<mat_t<double>> ProjectorAmplitudes(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block,
                                                           const mat_t<double>& C) const = 0;
    virtual std::vector<mat_t<dcmplx>> ProjectorAmplitudes(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block,
                                                           const mat_t<dcmplx>& C) const = 0;
    //! \brief Manifold \a M's U slots as projectors on its projector functions (the rows of \c ProjectorAmplitudes):
    //! a client splits a channel into the subspaces the +U term can give their own U.  THROWS before the term has
    //! seen a basis (the slots are built with its blocks).
    virtual HubbardSlots SlotProjectors(size_t M) const = 0;
};

} // namespace

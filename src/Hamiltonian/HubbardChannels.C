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
};

} // namespace

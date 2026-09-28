// File: Response/Probe.C  What a linear response PERTURBS and MEASURES, and the channel response matrix
// assembled from it (doc/LinearResponsePlan.md §2, stage R0).
//
// THE PAIR.  A channel J is perturbed with \f$\alpha\hat P_Je^{i\mathbf q\cdot\mathbf R}\f$ and every channel
// I is read as its occupation \f$n_I={\rm Tr}\hat P_I\rho\f$.  For a projector channel both halves are made
// of the SAME amplitudes \f$\ell=\langle w|\psi\rangle\f$: the perturbation's matrix on the block pair
// (k, k+q) is \f$A^J=\ell^{J\dagger}_{k+q}\ell^J_k\f$ and the measurement pairs \f$A^I\f$ with the first-order
// density matrix.  So \f$\chi_{IJ}=\langle A^I,\mathcal R\,A^J\rangle\f$ -- the Adjoint and Forward halves of
// one projector, which is what makes projector consistency hold by construction (insight 4).
//
// THE CONVENTION, stated once because the off-diagonal elements depend on it.  The perturbation is summed
// over LATTICE translations only, \f$\sum_{\mathbf R}e^{i\mathbf q\cdot\mathbf R}\hat P_{J,\mathbf R}\f$ -- the
// Bloch gauge of the basis (\f$e^{i\mathbf k\cdot\mathbf R_n}\f$).  A code that phases by \f$\mathbf R+\tau_J\f$
// (atom positions) differs in \f$\chi_{IJ}(\mathbf q)\f$ by \f$e^{i\mathbf q\cdot(\tau_I-\tau_J)}\f$; the
// diagonal and every REAL-SPACE element \f$\chi_{I\mathbf R,J\mathbf R'}\f$ are convention-free, so those are
// what an oracle comparison uses.
module;
#include <iosfwd>
#include <memory>
#include <string>
#include <vector>
export module qchem.Response.Probe;
export import qchem.Response.Reference;

export namespace qchem::Response
{

//! \brief The perturb/measure pair of a set of CHANNELS.  One face whatever the channel is: a Hubbard
//! manifold today; a dipole, an irrep slot or a spin-antisymmetric (J) channel later (§4, §6 D2).
class ChannelProbe
{
public:
    virtual ~ChannelProbe() = default;
    virtual size_t      NumChannels() const = 0;
    virtual std::string Label(size_t I) const = 0;
    //! Channel \a J's perturbation, modulated by \a rule (a lattice \f$e^{i\mathbf q\cdot\mathbf R}\f$, or
    //! \c Invariant), as a first-order Fock change on \a rule's block pairs (orbital basis).
    virtual BlockPairs Perturbation(size_t J, const SelectionRule& rule) const = 0;
    //! Every channel's first-order response \f$\delta n_I\f$ (per cell) from a first-order density on \a rule's pairs.
    virtual cvec_t     Measure(const SelectionRule& rule, const BlockPairs& dD) const = 0;
};

//! \brief Channels given by the orbitals' AMPLITUDES on each channel's projector functions,
//! \f$\hat P_J=\sum_\mu|w_\mu\rangle\langle w_\mu|\f$: \c amp[b][J] is \f$m_J\times n_{\rm orb}(b)\f$, the
//! \f$\ell=\langle w_\mu|\psi_i\rangle\f$ of block \a b's orbitals, in the reference's block order.
class AmplitudeProbe : public ChannelProbe
{
public:
    AmplitudeProbe(const Reference& ref, std::vector<std::vector<cmat_t>> amp, std::vector<std::string> labels);
    virtual size_t      NumChannels() const override {return itsLabels.size();}
    virtual std::string Label(size_t I) const override {return itsLabels[I];}
    virtual BlockPairs  Perturbation(size_t J, const SelectionRule& rule) const override;
    virtual cvec_t      Measure(const SelectionRule& rule, const BlockPairs& dD) const override;
private:
    const Reference&                 itsRef;
    std::vector<std::vector<cmat_t>> itsAmp;
    std::vector<std::string>         itsLabels;
};

//! \brief Channels given by one-body OPERATORS \f$\hat O_J\f$ (a dipole component, later an SOC or field
//! operator): channel J perturbs with \f$\hat O_J\f$ and every channel I reads \f$\langle\hat O_I\rangle\f$ --
//! the Adjoint/Forward pair again, with the SAME matrices both ways.  \c ops[J] is \f$\hat O_J\f$ in the
//! orbital basis on \a rule's block pairs (the \c OrbitalFrame's ToMO of its AO matrices).
class OperatorProbe : public ChannelProbe
{
public:
    OperatorProbe(const Reference& ref, std::vector<BlockPairs> ops, std::vector<std::string> labels,
                  std::shared_ptr<const SelectionRule> rule);
    virtual size_t      NumChannels() const override {return itsLabels.size();}
    virtual std::string Label(size_t I) const override {return itsLabels[I];}
    //! THROWS if \a rule pairs the blocks differently from the rule the operators were built on.
    virtual BlockPairs  Perturbation(size_t J, const SelectionRule& rule) const override;
    virtual cvec_t      Measure(const SelectionRule& rule, const BlockPairs& dD) const override;
private:
    void CheckRule(const SelectionRule& rule) const;
    const Reference&                     itsRef;
    std::vector<BlockPairs>              itsOps;
    std::vector<std::string>             itsLabels;
    std::shared_ptr<const SelectionRule> itsRule;
};

//! \brief The channel response matrix over a q-mesh: \f$\chi_{IJ}(\mathbf q)=\delta n_I/\delta\alpha_J\f$
//! (Hartree\f$^{-1}\f$), with the response gap it was gated on.
struct ChannelResponse
{
    std::vector<std::string> labels;
    ivec3_t                  Nq;
    std::vector<MeshShift>   q;       //!< the q-mesh, in CommensurateShifts order
    std::vector<cmat_t>      chi;     //!< per q, \f$n_{\rm ch}\times n_{\rm ch}\f$
    double                   gap=0;   //!< the smallest response gap over every q (inf when ungated)
    double                   noise=0; //!< the reference's eigenvalue noise that gap was gated against

    //! \brief The REAL-SPACE response over the q-mesh's supercell: row/column \f$(I,\mathbf R)\f$ at index
    //! \f$r\,n_{\rm ch}+I\f$ for image \f$\mathbf R=(r_x,r_y,r_z)\in[0,N_q)\f$ (r = (r_x N_y + r_y)N_z + r_z),
    //! \f$\chi_{I\mathbf R,J\mathbf R'}=N_q^{-1}\sum_{\mathbf q}e^{i\mathbf q\cdot(\mathbf R-\mathbf R')}\chi_{IJ}(\mathbf q)\f$ --
    //! the response of an ISOLATED perturbation, exactly what a supercell of that size would give (§1).
    //! THROWS if the result is not real to 1e-10 relative (time reversal is broken: a defect, not physics
    //! here -- no spin-orbit, no field).
    rmat_t RealSpace() const;
    //! The one-line-per-q report plus the real-space on-site block, in units of \a perUnit (27.2114 for
    //! eV\f$^{-1}\f$, hp.x's unit), and the gap.
    std::ostream& Write(std::ostream&, double perUnit=1.0, const std::string& unitName="1/Ha") const;
};

//! \brief \f$\chi_0\f$: the independent-particle channel response over the \a Nq q-mesh -- R0 alone, no
//! kernel (stage R0).  FAILS on an incommensurate q-mesh or an inverted/unresolved coupled pair (E1).
Outcome<ChannelResponse,ResponseFailure> IndependentResponse(const Reference&, const ChannelProbe&, ivec3_t Nq);

} // namespace

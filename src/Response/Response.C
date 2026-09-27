// File: Response/Response.C  qchem.Response -- the linear-response library's front door, and the ADAPTERS
// from a converged Bloch wave function to the library's own vocabulary (doc/LinearResponsePlan.md, stage R0).
//
// The Reference and the probes are built from VALUES (ReferenceBlock, amplitude matrices) so a unit test can
// feed them a model; these adapters are the one place that walks a real wave function to fill them.  They
// read only the WaveFunction's const face -- documented there as the face for "future property/post-HF
// code" -- and the Hamiltonian's HubbardChannels capability: nothing below qcWaveFunction learns about
// linear response.
module;
#include <memory>
#include <vector>
export module qchem.Response;
export import qchem.Response.Reference;
export import qchem.Response.Probe;
export import qchem.WaveFunction;                       // cWaveFunction
export import qchem.Hamiltonian.HubbardChannels;        // the +U projectors

export namespace qchem::Response
{

//! Which blocks share one chemical potential -- the reservoir partition the ground state was filled with
//! (Crystal_EC: \c globalFermi shares across k, \c spinsShareFermi across the two spin channels).  Only a
//! smeared rule's q = 0 Fermi shift reads it.
struct Reservoirs
{
    bool acrossK    = false;
    bool acrossSpin = false;
};

//! \brief The unperturbed state of a converged Bloch run: every block's eigenvalues and fractional
//! occupations (virtuals included), its BZ weight, over the run's own occupancy rule.
//! \a occ is the configuration the FINAL SCF stage filled with (the rule is rebuilt from the same value,
//! E1); \a eigenNoise is that stage's measured eigenvalue noise (Hartree).  THROWS on an IBZ-reduced run (D5).
Reference MakeReference(const WaveFunction::cWaveFunction& wf, const OccupationConfig& occ, Reservoirs res,
                        double eigenNoise);

//! \brief The run's Hubbard manifolds as probe channels: each block's orbitals' amplitudes on each channel's
//! projector functions, from the +U term's OWN projectors (\a hub, \c tHamiltonian::GetHubbardChannels) --
//! the same projector the ground state's occupations were made with.  \a ref must be \a wf's reference;
//! the probe keeps a reference to it.
AmplitudeProbe MakeHubbardProbe(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                                const Hamiltonian::HubbardChannels& hub);

} // namespace

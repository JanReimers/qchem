// File: Hamiltonian/Internal/ACBN0.C  ACBN0 -- (U-bar, J-bar) per Hubbard manifold from OUR OWN bare on-site
// two-electron integrals and a RENORMALISED density matrix (Agapito, Curtarolo, Buongiorno Nardelli, PRX 5,
// 011006 (2015), arXiv:1406.3259; doc/OpenWork.md step 5 increment 3; pin 23).
//
// THE IDEA.  Anisimov's on-site Hartree-Fock energy of the manifold (eq 1) is evaluated with BARE integrals
// (m1 m2|m3 m4) over the manifold's functions in the central cell and a density matrix in which every
// orbital is weighted by its own charge ON the manifold set: the electrons that sit elsewhere do not
// repel there.  That weighting is the screening; no Vee is ever chosen.  Comparing with Dudarev's averaged
// form (eq 2) defines
//     U-bar = Sum_{1234} (Pa+Pb)_{12}(Pa+Pb)_{34} (12|34)  /  [ (Na+Nb)^2 - Sum_m (Na_m^2 + Nb_m^2) ]     (eq 12)
//     J-bar = Sum_{1234} [Pa_{12}Pa_{34} + Pb_{12}Pb_{34}] (14|32)  /  [ Na^2 - Sum_m Na_m^2 + (b) ]      (eq 13)
// with P^sigma the RENORMALISED density matrix of the manifold (eq 10a, each orbital weighted by N-bar_i) but
// N^sigma_m the UNRENORMALISED populations (eq 10c carries no N-bar) -- the asymmetry that screens: U-bar
// scales as N-bar^2 for a manifold the KS states only partly live in.  U_eff = U-bar - J-bar is what the
// +U term takes.
//
// TWO BASES, AND THEY MUST NOT BE MIXED (found 2026-09-21: pairing Löwdin-basis coefficients with AO-basis
// integrals gave U-bar = 182 eV on MnO).  The integrals are over the RAW AOs phi_m of the manifold, so the
// density matrix in the numerator is the AO-basis one restricted to the manifold's functions -- the paper's
// eq 9, P-bar^sigma_{mm'} = Sum_k w_k Sum_i f_ki N-bar_ki c_mi c^*_m'i with the RAW coefficients c_mi of the
// manifold columns: its energy is the self-Coulomb of the d-AO component of the density.  LÖWDIN, NOT
// MULLIKEN, enters in the two CHARGES (user, 2026-09-16: Mulliken charges are basis-sensitive on a diffuse
// span -- the whole 136-span story): the RENORMALISED occupation of orbital i, N-bar_i = Sum over the
// manifolds of the SAME (species, l) of |T^dagger c_i|^2 (T = S^{1/2}[:,M], the +U projector's own; the
// paper's {m-bar}: both Mn on AFM-II MnO), and the per-function populations N^sigma_m in the pair-count
// denominators, the diagonal of the Löwdin matrix n-bar = Sum w f N-bar ell ell^dagger (the paper uses the
// Mulliken (PS)_mm there, eq 10c).  With N-bar_i == 1 the Löwdin matrix is exactly the +U occupation matrix
// n -- the estimator also accumulates the unweighted sums and reports the bare averages, so the screening
// the renormalisation supplied is visible.
//
// WHAT IT IS NOT (yet): the true variational ACBN0 functional (U depending on the density inside the SCF);
// the paper's practice is the OUTER LOOP -- SCF at U^(n), estimate U^(n+1), repeat to 1e-4 eV -- and that
// is what the facade drives.  Orbital resolution (a U-bar per site-irrep slot) is a diagnostic to add on
// top: the same sums restricted to a slot's eigenvectors.
module;
#include <iosfwd>
#include <map>
#include <memory>
#include <vector>
export module qchem.Hamiltonian.Internal.ACBN0;
export import qchem.Hamiltonian.HubbardEstimator;   // the face this realises
export import qchem.Hamiltonian.Internal.Hubbard;   // HubbardProjection (the term face the composite casts to), HubbardManifold
import qchem.BasisSet.BareCoulombSource;            // ERI4Block
import qchem.Symmetry.Irrep;                        // Spin
import qchem.Blaze;

export namespace qchem::Hamiltonian
{

class ACBN0 : public virtual HubbardUEstimator
{
public:
    //! \a term outlives the estimator (the Hamiltonian owns it; the facade owns the Hamiltonian).
    explicit ACBN0(HubbardProjection& term);
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const Spin& s, double w,
                            const mat_t<double>& C, const rvec_t& f) override;
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const Spin& s, double w,
                            const mat_t<dcmplx>& C, const rvec_t& f) override;
    virtual std::vector<HubbardEstimate> Evaluate() const override;
    virtual std::ostream& Write(std::ostream&) const override;
    virtual void Apply(const std::vector<HubbardEstimate>& e) override;

private:
    template <class U> void AccumulateT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block, const Spin& s, double w,
                                        const mat_t<U>& C, const rvec_t& f);
    //! The bare integrals of manifold \a M, from the first block that carries the face (geometry-fixed).
    template <class U> void EnsureIntegrals(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block);
    //! The two averages (eqs 12, 13): AO-basis matrices \a Pa/\a Pb in the numerators, the populations
    //! \a Na/\a Nb in the pair-count denominators.
    static void Averages(const rmat_t& Pa, const rmat_t& Pb, const rvec_t& Na, const rvec_t& Nb,
                         const BasisSet::ERI4Block& eri, double& Ubar, double& Jbar);

    HubbardProjection&                     itsTerm;
    std::vector<BasisSet::ERI4Block>       itsERI;       //!< per manifold (empty until the first block)
    //! Per manifold: the AO-basis density matrix on the manifold's functions (the numerator) and the Löwdin
    //! occupation matrix (its diagonal = the populations of the denominators), each renormalised and unweighted.
    struct Channel { std::vector<rmat_t> P, Pbare, L, Lbare; };
    std::map<Spin,Channel>                 itsChannels;  //!< Up/Down (a Spin::None block feeds both, halved)
    bool                                   itsFed=false;
};

} // namespace

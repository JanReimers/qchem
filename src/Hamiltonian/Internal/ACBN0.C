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
// with P^sigma the renormalised density matrix of the manifold, N^sigma_m = P^sigma_mm, N^sigma = Tr P^sigma,
// and U_eff = U-bar - J-bar what the +U term takes.
//
// LÖWDIN, NOT MULLIKEN (user, 2026-09-16: Mulliken charges are basis-sensitive on a diffuse span -- the
// whole 136-span story).  Per orbital i on block k: the Löwdin coefficients ell_i = T^dagger c_i in every
// manifold (T = S^{1/2}[:,M], the +U projector's own), the RENORMALISED occupation N-bar_i = Sum over the
// manifolds of the SAME (species, l) of |ell_i|^2 (the paper's {m-bar}: both Mn on AFM-II MnO), and
//     P-bar^sigma_M = Sum_k w_k Sum_i f_ki N-bar_ki ell_ki ell_ki^dagger .
// With N-bar_i == 1 this is exactly the +U occupation matrix n -- the estimator also accumulates that
// unweighted sum and reports the bare averages, so the screening the renormalisation supplied is visible.
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
    explicit ACBN0(const HubbardProjection& term);
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const Spin& s, double w,
                            const mat_t<double>& C, const rvec_t& f) override;
    virtual void Accumulate(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const Spin& s, double w,
                            const mat_t<dcmplx>& C, const rvec_t& f) override;
    virtual std::vector<HubbardEstimate> Evaluate() const override;
    virtual std::ostream& Write(std::ostream&) const override;

private:
    template <class U> void AccumulateT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block, const Spin& s, double w,
                                        const mat_t<U>& C, const rvec_t& f);
    //! The bare integrals of manifold \a M, from the first block that carries the face (geometry-fixed).
    template <class U> void EnsureIntegrals(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& block);
    //! The two averages (eqs 12, 13) of a pair of channel matrices.
    static void Averages(const rmat_t& Pa, const rmat_t& Pb, const BasisSet::ERI4Block& eri, double& Ubar, double& Jbar);

    const HubbardProjection&               itsTerm;
    std::vector<BasisSet::ERI4Block>       itsERI;       //!< per manifold (empty until the first block)
    struct Channel { std::vector<rmat_t> P, Pbare; };   //!< per manifold: renormalised and unweighted
    std::map<Spin,Channel>                 itsChannels;  //!< Up/Down (a Spin::None block feeds both, halved)
    bool                                   itsFed=false;
};

} // namespace

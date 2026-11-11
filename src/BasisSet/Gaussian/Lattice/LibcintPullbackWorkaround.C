// File: BasisSet/Gaussian/Lattice/LibcintPullbackWorkaround.C  ⚠ TEMPORARY: the GPW (MnD) lattice-sum overlap and kinetic DERIVATIVES from
// libcint's ip integrals, while Enzyme cannot yet differentiate the MnD kernel (EnzymeAD/Enzyme#3747, a silently wrong reverse-mode gradient,
// whose trigger -- a callee filling a loop-local array through an inner loop with a data-dependent trip count -- is the MnD Hermite tables).
//
// THE USER'S RULING (2026-11-10): "libCint and our own MnD code should not be entangled (outside oracle unit tests).  We should be able to
// drop one without affecting the other.  Having said all that I don't mind a temporary work around while waiting for the Enzyme fix."  So the
// entanglement lives HERE and nowhere else: PG_Cart::Orbital_IBS's two LatticeSum1E pullback bodies are its only callers (one line each), and
// they hand it MnD's OWN screened image enumerator (ForImageOffsets: the energy series' list, never a second screen -- pin 1).  When #3747 is
// fixed, those two bodies become Enzyme over the MnD kernel and this module is deleted (doc/OpenWork.md row "replace with Enzyme over MnD").
// The libcint per-image integrals it calls (PG_LibCint::NR_Evaluator::ShellPairImage*) are libcint's own faces and stay.
// The fix is underway upstream: EnzymeAD/Enzyme PR #3753 addresses #3747 (and #3748).  Once it is merged, ONE Enzyme pin move takes both
// fixes, and Enzyme over MnD replaces this seam.
//
// While it lives, two gates hold it (UTGaussian_BS M_LatticePullback): the gradient and virial against the central difference of MnD's
// ENERGY-path lattice sum, and libcint's value integral against MnD's image by image (so a mismatch is attributed to the right code).
module;
#include <array>
#include <complex>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.BasisSet.Gaussian.Lattice.LibcintPullbackWorkaround;
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;              // GeometryPullback, cellphase_t
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;     // the MnD block's component table (the libcint twin is built on it)
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF; // a component's centre
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;             // ShellPairImageGradient: <∂p|O|q(· − R)>
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

export namespace qchem::BasisSet::Gaussian::LibcintWorkaround
{

//! The energy series' screened image list for the ordered component pair (i, j): calls back (n, R_n) for every kept offset.
using ImageEnumerator = std::function<void(size_t i, size_t j, const std::function<void(const ivec3_t& n, const rvec3_t& Roff)>& cb)>;

//! ⚠ TEMPORARY pending EnzymeAD/Enzyme#3747 (the user's 2026-11-10 ruling, file head).  The gradient and image virial of
//! E = Re Tr(Ā† M), M = MakeOverlap(phase,A) (\a which 0) or MakeKinetic(phase,A) (\a which 1) of the MnD block \a pg (LatticeSum1E's doc).
//!
//! The series is M_ij = Σ_n phase(n) K_ij(n), K_ij(n) = <i|O|j(· − R_n)> with i on site a at R_i and j on site b at R_j.  Over the FULL
//! Hermitian matrix E = Σ_ij Re(Ā_ij* M_ij) (the series fills the upper triangle; the lower is its conjugate, the ordered pair (j,i) walking
//! −n).  libcint's X_c(i,j_n) = <∂_c i|O|j_n> is the electron-coordinate derivative of the bra, so ∂K/∂R_a = −X and, by translation
//! invariance, ∂K/∂(R_j+R_n) = +X.  With w_c = Re(Ā_ij* phase(n) X_c): g_a −= w, g_b += w, and the virial over the two unfolded centres is
//! −w ⊗ R_i + w ⊗ (R_j + R_n) = −w ⊗ d, d = R_i − R_j − R_n, which carries no origin (Σ g = 0 term by term).
inline LatticeSum1E::GeometryPullback LatticeOneElectronPullback(int which, const Evaluators::PG_Cart_MnD::PGData& pg, const ImageEnumerator& images,
                                                                 const chmat_t& Abar, const LatticeSum1E::cellphase_t& phase, const Structure* st)
{
    const char* what = which==0 ? "the lattice overlap pullback" : "the lattice kinetic pullback";
    if (!st) throw std::invalid_argument(std::string(what)+": needs the structure (the sites the functions sit on)");
    if (Abar.rows()!=pg.size()) throw std::invalid_argument(std::string(what)+": the adjoint is "+std::to_string(Abar.rows())+
                                                            " square for "+std::to_string(pg.size())+" functions");
    const Evaluators::PG_LibCint::NR_Evaluator lc(pg, st);
    LatticeSum1E::GeometryPullback r{rvec3vec_t(st->GetNumAtoms(), rvec3_t(0,0,0)), rmat3d_t(0,0,0, 0,0,0, 0,0,0)};
    double V[3][3]={{0,0,0},{0,0,0},{0,0,0}};
    for (size_t si=0;si<lc.NumShells();si++)
    {
        const std::vector<size_t>& P=lc.ShellComponents(si);
        const size_t a=lc.ShellAtom(si);
        const rvec3_t Ri=pg.radials[P[0]]->GetCenter();
        for (size_t sj=0;sj<lc.NumShells();sj++)
        {
            const std::vector<size_t>& Q=lc.ShellComponents(sj);
            const size_t b=lc.ShellAtom(sj);
            const rvec3_t Rj=pg.radials[Q[0]]->GetCenter();
            images(P[0], Q[0], [&](const ivec3_t& n, const rvec3_t& Roff)   // a shell pair shares its radials, so one list serves it
            {
                const std::array<rmat_t,3> X=lc.ShellPairImageGradient(which, si, sj, Roff);
                const dcmplx ph=phase(n);
                double w[3]={0,0,0};
                for (size_t p=0;p<P.size();p++)
                    for (size_t q=0;q<Q.size();q++)
                    {
                        const dcmplx c=std::conj(dcmplx(Abar(P[p],Q[q])))*ph;
                        for (int k=0;k<3;k++) w[k]+=std::real(c*X[k](p,q));
                    }
                const rvec3_t wv(w[0],w[1],w[2]), d=Ri-Rj-Roff;
                r.gradient[a]-=wv;
                r.gradient[b]+=wv;
                const double dd[3]={d.x,d.y,d.z};
                for (int u=0;u<3;u++) for (int v=0;v<3;v++) V[u][v]-=w[u]*dd[v];
            });
        }
    }
    r.virial=rmat3d_t(V[0][0],V[0][1],V[0][2], V[1][0],V[1][1],V[1][2], V[2][0],V[2][1],V[2][2]);
    return r;
}

} // namespace qchem::BasisSet::Gaussian::LibcintWorkaround

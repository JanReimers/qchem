// File: BasisSet/SeparablePotential.C  Separable (Kleinman-Bylander) NONLOCAL potential.
//
// Rung 2 of "lineage A" (see doc/PlaneWavePlan.md): a norm-conserving pseudopotential's nonlocal part
// in Kleinman-Bylander separable form,
//     V_NL = Sum_{atom a} Sum_{projector p} |beta^a_p> D_p <beta^a_p| .
// In a plane-wave basis each projector contributes a reciprocal-space radial form factor
// beta-tilde_p(|q|) (q = k+G) and a KB coefficient D_p; the basis set (PlaneWave_IBS::
// MakeSeparablePotential) folds in 1/Omega and the structure-factor phase e^{-i(G-G').tau_a}.
//
// This is the nonlocal sibling of LocalPotential: the external one-body potential is V = V_loc +
// V_nonlocal, both model-parameterized contributions to the SAME external block.
//
// Angular channels: each projector carries an angular momentum l (AngularMomentum); MakeSeparablePotential
// weights its rank-1 radial product by the addition-theorem factor (2l+1) P_l(cos gamma), gamma the angle
// between k+G and k+G' -- the SAME angular structure the APW/LAPW sphere terms use.  l=0 (P_0=1) is the
// spherically-symmetric s-channel.  The radial form factor beta-tilde_l(|q|) is a model input here (the
// production norm-conserving / PAW route obtains it from a spherical-Bessel transform j_l of beta_l(r)).
module;
#include <cassert>
#include <cmath>
#include <utility>
#include <vector>

export module qchem.BasisSet.SeparablePotential;
import qchem.Math;   // Pi, FourPi

export namespace BasisSet
{

//! \brief A separable Kleinman-Bylander nonlocal potential, supplying per-species projector form
//! factors and KB coefficients.  Open/closed extension point for the nonlocal part of a pseudopotential.
class SeparablePotential
{
public:
    virtual ~SeparablePotential() {}
    //! Number of (radial) projectors for nuclear species \a Z.
    virtual size_t NumProjectors(int Z) const=0;
    //! Kleinman-Bylander coefficient \f$D_p\f$ for projector \a p of species \a Z.  [energy]
    virtual double Coefficient  (int Z, size_t p) const=0;
    //! Reciprocal-space radial projector form factor \f$\tilde\beta_p(|q|)\f$, \f$q=|k+G|\f$.
    virtual double Projector    (int Z, size_t p, double q) const=0;
    //! Angular momentum \a l of projector \a p (default 0 = s-channel).  Sets the \f$(2l+1)P_l(\cos\gamma)\f$
    //! angular weight of this projector's contribution.
    virtual int    AngularMomentum(int /*Z*/, size_t /*p*/) const {return 0;}
};

//! \brief A single Gaussian Kleinman-Bylander projector in channel \a l,
//! \f$\tilde\beta(q)=e^{-\sigma^2 q^2/2}\f$ with coefficient \f$D\f$ -- a minimal analytic demonstrator
//! of the separable nonlocal structure (each atom contributes a rank-1 \f$V_{NL}\f$ per channel).
//! Defaults to l=0 (the s-channel), so existing callers are unaffected.
class GaussianProjector : public SeparablePotential
{
public:
    GaussianProjector(double sigma, double D, int l=0) : itsSigma(sigma), itsD(D), itsL(l) {}
    virtual size_t NumProjectors(int) const {return 1;}
    virtual double Coefficient  (int, size_t) const {return itsD;}
    virtual double Projector    (int, size_t, double q) const {return std::exp(-0.5*itsSigma*itsSigma*q*q);}
    virtual int    AngularMomentum(int, size_t) const {return itsL;}
private:
    double itsSigma; //!< Projector width (Bohr).
    double itsD;     //!< KB coefficient (energy).
    int    itsL;     //!< Angular-momentum channel.
};

//! \brief A real norm-conserving pseudopotential's nonlocal part in the analytic Goedecker / Hartwigsen-
//! Goedecker-Hutter (HGH) form.  Each angular channel \a l supplies a projector radius \f$r_l\f$ and a
//! symmetric coupling matrix \f$h^l_{ij}\f$ (1-3 projectors).  The reciprocal radial projectors are the
//! closed-form HGH momentum-space functions \f$\tilde\beta_i^l(q)=\frac{1}{4\pi}\,
//! \pi^{5/4} q^l \sqrt{r_l^{2l+3}}\,Q_i^l(qr_l)\,e^{-(qr_l)^2/2}\f$ (the \f$4\pi\f$ makes them the bare
//! spherical-Bessel transforms \f$\int\beta_i^l(r)j_l(qr)r^2dr\f$, cross-checked analytically).
//!
//! The coupled \f$h^l\f$ is diagonalised once into independent Kleinman-Bylander projectors: with
//! \f$h^l=\sum_\alpha\lambda_\alpha v_\alpha v_\alpha^T\f$, projector \f$\alpha\f$ has KB coefficient
//! \f$D_\alpha=\lambda_\alpha\f$ and form factor \f$\tilde\beta_\alpha(q)=\frac{1}{\sqrt{4\pi}}
//! \sum_i v_{\alpha,i}\,\pi^{5/4} q^l\sqrt{r_l^{2l+3}}Q_i^l e^{-(qr_l)^2/2}\f$, so the generic
//! (2l+1)P_l(cosγ) assembler in PlaneWave_IBS reproduces \f$\frac1\Omega(2l+1)P_l\sum_{ij}\tilde\beta_i
//! h_{ij}\tilde\beta_j\f$ exactly.
class HGH_SeparablePotential : public SeparablePotential
{
public:
    //! Build by adding channels (AddChannel) from real GTH parameters; the GTH database reader
    //! (GetGTH in GTH_Potentials.C) is the per-element source, replacing hardcoded factories.

    virtual size_t NumProjectors  (int)         const {return itsProj.size();}
    virtual double Coefficient    (int, size_t p) const {return itsProj[p].D;}
    virtual int    AngularMomentum(int, size_t p) const {return itsProj[p].l;}
    virtual double Projector      (int, size_t p, double q) const
    {
        const Proj& pr=itsProj[p];
        double s=0.0;                                          // (1/sqrt 4pi) Sum_i v_i projG_i(q)
        for (size_t i=0;i<pr.v.size();i++) s += pr.v[i]*ProjG(q, pr.l, static_cast<int>(i), pr.rl);
        return s/std::sqrt(FourPi);
    }

    //! \brief Add an angular channel from its symmetric KB coefficient matrix \a h (the HGH h-matrix
    //! for momentum \a l, with projector radius \a rl): diagonalise into independent KB projectors
    //! \f$|\beta\rangle D\langle\beta|\f$.  A 0x0 \a h (a tabulated channel with no projectors) adds
    //! nothing.  Public so the GTH table reader can assemble a potential channel-by-channel.
    void AddChannel(int l, double rl, const std::vector<std::vector<double>>& h)
    {
        for (auto& [lambda,v] : SymEig(h)) itsProj.push_back(Proj{l, rl, lambda, std::move(v)});
    }

private:
    struct Proj { int l; double rl; double D; std::vector<double> v; };  //!< one diagonalised KB projector
    std::vector<Proj> itsProj;

    //! Eigenpairs (lambda, normalised eigenvector) of a small symmetric matrix, by cyclic Jacobi.  The KB
    //! form needs the spectral decomposition Sum_p lambda_p v_p v_p^T, which is invariant to eigenvector
    //! sign and ordering -- so any correct symmetric eigensolver gives the same potential.  Handles n=0
    //! (no projectors) and n=1 (trivial) as well as the 2x2/3x3 HGH channels.
    static std::vector<std::pair<double,std::vector<double>>> SymEig(const std::vector<std::vector<double>>& h)
    {
        size_t n=h.size();
        std::vector<std::pair<double,std::vector<double>>> out;
        if (n==0) return out;

        std::vector<std::vector<double>> A=h;                              // working copy, diagonalised in place
        std::vector<std::vector<double>> V(n, std::vector<double>(n,0.0)); // accumulated eigenvectors (-> identity)
        for (size_t i=0;i<n;i++) V[i][i]=1.0;

        for (int sweep=0; sweep<100; sweep++)
        {
            double off=0.0;
            for (size_t p=0;p<n;p++) for (size_t q=p+1;q<n;q++) off+=A[p][q]*A[p][q];
            if (off<1e-30) break;
            for (size_t p=0;p<n;p++) for (size_t q=p+1;q<n;q++)
            {
                if (std::abs(A[p][q])<1e-300) continue;
                double theta=(A[q][q]-A[p][p])/(2*A[p][q]);
                double t=(theta>=0?1.0:-1.0)/(std::abs(theta)+std::sqrt(theta*theta+1.0));
                double c=1.0/std::sqrt(t*t+1.0), s=t*c;
                for (size_t k=0;k<n;k++) { double a=A[k][p],b=A[k][q]; A[k][p]=c*a-s*b; A[k][q]=s*a+c*b; }
                for (size_t k=0;k<n;k++) { double a=A[p][k],b=A[q][k]; A[p][k]=c*a-s*b; A[q][k]=s*a+c*b; }
                for (size_t k=0;k<n;k++) { double a=V[k][p],b=V[k][q]; V[k][p]=c*a-s*b; V[k][q]=s*a+c*b; }
            }
        }
        for (size_t i=0;i<n;i++)
        {
            std::vector<double> vec(n);
            for (size_t k=0;k<n;k++) vec[k]=V[k][i];
            out.push_back({A[i][i], std::move(vec)});
        }
        return out;
    }

    //! HGH momentum-space polynomial Q_i^l(x), x=q r_l (Goedecker/HGH; pyscf _qli, 0-indexed i).
    static double Qli(double x, int l, int i)
    {
        using std::sqrt;
        if (l==0 && i==0) return 4*sqrt(2.0);
        if (l==0 && i==1) return 8*sqrt(2.0/15)*(3-x*x);
        if (l==0 && i==2) return (16.0/3)*sqrt(2.0/105)*(15-10*x*x+x*x*x*x);
        if (l==1 && i==0) return 8*sqrt(1.0/3);
        if (l==1 && i==1) return 16*sqrt(1.0/105)*(5-x*x);
        if (l==2 && i==0) return 8*sqrt(2.0/15);
        assert(false && "HGH Q_i^l: channel not tabulated"); return 0.0;
    }

    //! Reciprocal HGH projector projG_li(q) = Q_i^l(q r_l) pi^{5/4} q^l sqrt(r_l^{2l+3}) e^{-(q r_l)^2/2}.
    static double ProjG(double q, int l, int i, double rl)
    {
        double x=q*rl, ql=(l==0)?1.0:std::pow(q, static_cast<double>(l));
        return Qli(x,l,i)*std::pow(Pi,1.25)*ql*std::sqrt(std::pow(rl,2*l+3))*std::exp(-0.5*x*x);
    }
};

} //namespace

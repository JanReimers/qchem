// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/Evaluator.C
//
// First concrete molecular evaluator (Goal B): the NON-RELATIVISTIC, CARTESIAN polarized-Gaussian basis
// integrated by McMurchie-Davidson (hence PG_Cart_MnD::NR_Evaluator).  IS-A PolarizedGaussian PGData --
// the flattened (radial x polarization) component set -- and is meant to be a base subobject of the PG
// IrrepBasisSet (the atom evaluators relate to their IBS the same way).  The inline 1E kernels are thin
// wrappers over GaussianRF::Integrate, the very primitives that M_PG_Oracle validates to machine
// precision, with the per-component normalization n_i*n_j folded in as PolarizedGaussian's MakeIntegrals.
//
// Satisfies is1E_Evaluator, isDFT_Evaluator and isHF_Evaluator (checked below).  The VectorFunction
// grid-eval (operator()(r), needed before isOpr can fold into isDFT) comes in a later increment.
//
// Sibling evaluators anticipated: PG_Spherical_MnD, PG_Cart_PRISM, ... each a new Evaluator over (mostly)
// the same shared PGData -- which is why PGData is kept separate, not absorbed here.
module;
#include <cassert>
#include <complex>   // std::real (Hermitian-diagonal projection of the Bloch lattice sum)
#include <string>
#include <ostream>
#include <iostream>   // std::cerr (the one-line stream-cache coverage readout)
#include <vector>
#include <cmath>      // std::sqrt/std::log (the lattice-sum magnitude-screening reach radius)
#include <algorithm>  // std::min (alpha_min per radial)
#include <functional> // cellphase_t (the caller-supplied Bloch phase of a cross-cell offset)
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;      // PGData
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;  // GaussianRF named kernels
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;// Polarization
import qchem.BasisSet.Molecule.LatticeSum1E;                       // GaussianFunction (the <chi_i|g> seam type)
import qchem.Structure;
import qchem.UnitCell;                                             // UnitCell (ToCartesian/ToFractional: grid<->cell for collocation)
import qchem.Types;
import qchem.Blaze;                                                // rsmat_t (the lattice-sum matrices)

export namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

class NR_Evaluator : public virtual Evaluator, public PGData
{
public:
    // The evaluator OWNS its data (IS-A PGData) rather than viewing one: it is a base subobject of the PG
    // IrrepBasisSet, exactly as the atom evaluators are bases of their IBS.  The default constructor builds
    // an empty PGData; the host IBS fills it via PGData::Init.  The structure is NOT evaluator state -- it is
    // passed per call to the Nuclear kernel.
    NR_Evaluator() = default;

    // --- cold-path Evaluator interface ---
    virtual size_t        size    () const {return PGData::size();}
    virtual rvec_t        Norm    () const {return ns;}
    virtual std::string   Name    () const {return "PolarizedGaussian";}
    virtual std::ostream& Write   (std::ostream& os) const {return os << "PG_Cart_MnD::NR_Evaluator[" << size() << "]";}

    // --- 1E inline kernels (hot loops) ---
    // Grad2 is the FULL Cartesian \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block: NO 1/2 (that is
    // the Hamiltonian's), NO centrifugal term (atom-only).  Nuclear is the multi-centre attraction.
    double Norm   (size_t i)          const {return ns[i];}
    double Overlap(size_t i,size_t j) const {return radials[i]->Overlap2C(*radials[j], pols[i], pols[j]) * ns[i]*ns[j];}
    double Grad2  (size_t i,size_t j) const {return radials[i]->Grad2    (*radials[j], pols[i], pols[j]) * ns[i]*ns[j];}
    double Nuclear(size_t i,size_t j,const Structure* cl=0) const
    {
        return radials[i]->Nuclear(*radials[j], pols[i], pols[j], cl) * ns[i]*ns[j];
    }

    // --- Periodic (lattice-summed) 1E matrices at general k -- the GPW seam ------------------------------
    // Sum_R e^{ik.R} <chi_i | O | chi_j(.-R)>, Hermitian (real at Gamma where every phase is 1).  Rs are the
    // Cartesian lattice translations (MUST include {0}) and phases[r]=e^{ik.Rs[r]} the matching per-image
    // Bloch weight (origin phase 1); chi_j(.-R) is radials[j] placed at its centre + R via GaussianRF::
    // AtCenter.  The R=0 term is the finite kernel above, so Rs={0}/phase 1 reproduces Overlap/Grad2/Nuclear
    // exactly -- the SAME analytic M&D kernels, only the second centre shifted per image.  For an inversion-
    // symmetric Rs the sum is Hermitian, so we fill the upper triangle (the (i,j) Bloch element) and let
    // chmat_t mirror the lower as its conjugate; the diagonal is real by that same symmetry.  These realise
    // Molecule::LatticeSum1E (on the host IBS).  (Rs,phases) is a {R}+{e^{ik.R}} weighted point set -- FUTURE:
    // one const cMesh& (Mesh<dcmplx>) once qcMesh grows complex weights (see LatticeSum1E.C header).
    // MAGNITUDE SCREENING (CP2K's EPS_PGF_ORB) replaces a fixed geometric |R|<=Rcut cutoff.  Per component the
    // reach r_i = sqrt(-ln(eps)/alpha_min_i) is where its MOST DIFFUSE primitive falls below eps; a pair
    // (i home, j at c_j+R) whose centres are farther than r_i+r_j apart has a Gaussian PRODUCT -- hence every
    // 1E integrand (overlap, <p^2>, and the nuclear product.(1/r)) -- below eps, so its lattice term is dropped.
    // This (a) keeps the sum SPARSE (a generous enumeration Rcut costs ~nothing -- tight functions reach nothing)
    // and (b) keeps S POSITIVE-DEFINITE at any enumerated reach: it drops ONLY sub-eps terms, so ||S-S_exact||<eps
    // (<< lambda_min) and S_exact (the full Bloch Gram) is PSD -- no Gibbs from a sharp cutoff chopping a still-
    // significant diffuse tail (the old Rcut=2a / SR-basis crutch).  reach_i+reach_j is CONSERVATIVE (>= the exact
    // pair threshold sqrt(-ln eps.(1/ai+1/aj))), so no significant term is ever dropped.  eps=1e-10 (numerically
    // exact for GPW tolerances).  NOTE: the enumerated Rs must still REACH far enough (screening only removes, it
    // cannot add a far term the caller never enumerated) -- so pass a generous Rcut and let the screen prune it.
    static constexpr double kScreenEps=1e-10;
    //! D-AWARE density-magnitude screen (CP2K's eps/|coef| radii, doc/GPWPlan.md 0a).  What lands on the grid
    //! is weight*chi_i*chi_j (weight ~ the density-matrix element), and the accuracy target is ABSOLUTE on the
    //! density -- so the tolerance a pair's box must honour is kDensityEps/|weight|, not kScreenEps: a
    //! small-|D| pair keeps a SMALLER box (radius shrinks with sqrt(ln), work with its cube), and a
    //! (pair, offset) whose |weight|*max|value| is below kDensityEps is skipped whole.  Dropping only
    //! sub-eps DENSITY contributions with smooth exponential tails: a magnitude screen, no Gibbs (the
    //! ringing ledger's class).  Clamped so a |weight|>1 never grows a box beyond the geometry screen
    //! (replay is capped by what was built).  The same criterion drives the integrate-back when the caller
    //! supplies its density (the SHARED ACTIVE SET: both directions skip identical terms, so the
    //! collocate/integrate adjoint stays machine-exact on the kept operator).
    static constexpr double kDensityEps=1e-10;
    //! The Bloch phase of an integer cell offset (== Molecule::LatticeSum1E::cellphase_t -- the same
    //! std::function type; the k-CONVENTION stays with the lattice-side caller, Gamma = the constant 1).
    using cellphase_t = std::function<dcmplx(const ivec3_t& n)>;
    template <class Kernel> chmat_t LatticeSum(const std::vector<rvec3_t>& Rs, const cvec_t& phases, Kernel K) const
    {
        assert(phases.size()==Rs.size());
        const double lne=-std::log(kScreenEps);
        std::vector<double> reach(size());
        for (auto i:indices())
        {
            double amin=radials[i]->GetExponents()[0];
            for (double e:radials[i]->GetExponents()) amin=std::min(amin,e);
            reach[i]=std::sqrt(lne/amin);
        }
        chmat_t S(size());
        for (auto i:indices()) for (auto j:indices(i))
        {
            dcmplx s(0.0);
            const rvec3_t ci=radials[i]->GetCenter(); const rvec3_t cj=radials[j]->GetCenter();
            const double  rr=(reach[i]+reach[j])*(reach[i]+reach[j]);   // screen on squared distance (no sqrt)
            for (size_t r=0; r<Rs.size(); r++)
            {
                const rvec3_t d=ci-(cj+Rs[r]);
                if (d.x*d.x+d.y*d.y+d.z*d.z > rr) continue;             // magnitude screen: product < eps -> skip
                s += phases[r] * K(i, j, radials[j]->AtCenter(cj+Rs[r]));
            }
            s *= ns[i]*ns[j];
            S(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) auto-set to conj
        }
        return S;
    }
    chmat_t MakeOverlap(const std::vector<rvec3_t>& Rs, const cvec_t& phases) const
    {   return LatticeSum(Rs,phases,[this](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Overlap2C(cj,pols[i],pols[j]);}); }
    chmat_t MakeKinetic(const std::vector<rvec3_t>& Rs, const cvec_t& phases) const
    {   return LatticeSum(Rs,phases,[this](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Grad2   (cj,pols[i],pols[j]);}); }
    chmat_t MakeNuclear(const std::vector<rvec3_t>& Rs, const cvec_t& phases, const Structure* cl) const
    {   return LatticeSum(Rs,phases,[this,cl](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Nuclear(cj,pols[i],pols[j],cl);}); }

    // Molecule::LatticeSum1E: the lattice-summed overlap of every basis function with ONE caller-supplied
    // Cartesian-Gaussian function g -- b_i = Sum_R phases[R] <chi_i | g(.-Rs[R])>.  g stands in the chi_j
    // slot of the 2-centre kernel: each monomial term is a Polarization on an uncontracted GaussianRF, so
    // the sum reuses the analytic M&D Overlap2C verbatim.  Same magnitude screen as the pair sums (g's
    // reach from its own exponent); chi_i carries ns[i], g is integrated RAW (its scale is in g.terms).
    cvec_t MakeOverlap(const std::vector<rvec3_t>& Rs, const cvec_t& phases,
                       const Molecule::LatticeSum1E::GaussianFunction& g) const
    {
        assert(phases.size()==Rs.size());
        const double lne=-std::log(kScreenEps);
        const double reachG=std::sqrt(lne/g.alpha);
        int Lg=0;                                          // g's max polynomial degree (Hermite table sizing)
        for (const auto& t : g.terms) Lg=std::max(Lg, t.p.n+t.p.l+t.p.m);
        cvec_t b(size(), dcmplx(0.0));
        for (auto i:indices())
        {
            double amin=radials[i]->GetExponents()[0];
            for (double e:radials[i]->GetExponents()) amin=std::min(amin,e);
            const double reach=std::sqrt(lne/amin)+reachG;
            const rvec3_t ci=radials[i]->GetCenter();
            dcmplx s(0.0);
            for (size_t r=0; r<Rs.size(); r++)
            {
                const rvec3_t c=g.center+Rs[r];
                const rvec3_t d=ci-c;
                if (d.x*d.x+d.y*d.y+d.z*d.z > reach*reach) continue;   // product < eps -> skip
                GaussianRF gr(g.alpha, c, Lg);
                double sr=0.0;
                for (const auto& t : g.terms)
                    sr += t.c * radials[i]->Overlap2C(gr, pols[i], Polarization(t.p));
                s += phases[r]*sr;
            }
            b[i]=ns[i]*s;
        }
        return b;
    }

    // Molecule::LatticeSum1E: the finest exponent over every component's radial primitives -- the density
    // grid's resolution driver (see LatticeSum1E::MaxExponent).  Walks the owned radials; no primitive escapes.
    double MaxExponent() const
    {
        double amax=0.0;
        for (auto i:indices()) for (double e : radials[i]->GetExponents()) if (e>amax) amax=e;
        return amax;
    }

    // Molecule::LatticeSum1E: the COARSEST primitive exponent -- the diffuse end of the basis, mirroring
    // MaxExponent.  It sets the coarsest useful GPW density grid LEVEL (a product of two alpha_min primitives
    // has exponent 2*alpha_min, resolved by a proportionally coarse cutoff), so the level ladder runs from the
    // fine grid down to ~ this.  A scalar summary; no primitive escapes.
    double MinExponent() const
    {
        double amin=radials[0]->GetExponents()[0];
        for (auto i:indices()) for (double e : radials[i]->GetExponents()) if (e<amin) amin=e;
        return amin;
    }

    // --- ANALYTIC density collocation + integrate-back (CP2K GPW) -------------------------------------------
    // The periodic density is a product of BLOCH orbitals: rho = Sum_ij D_ij chi_i^k conj(chi_j^k), and
    // chi_i^k conj(chi_j^k) = Sum_R'' e^{-ik.R''} [ chi_i^0 . chi_j^R'' ] tiled -- so collocation loops the
    // CROSS-CELL pairs (chi_i home, chi_j in cell R'', phase from the caller), each collocated ANALYTICALLY on
    // its compact exp-tail box and MODULO-WRAPPED onto the grid.  TWO mechanisms, both screened, NO hard cutoff:
    //  - the cross-cell offset R (ForImageOffsets): MAGNITUDE-screened (include (i,j,R) only where chi_i and
    //    chi_j^R still overlap, reach_i+reach_j from the diffuse ends) -- like the 1E lattice sums, NOT a fixed
    //    Rcut, so no Gibbs ringing (screening drops only sub-eps pairs);
    //  - the modulo wrap (ForPairBox): tiles each compact box onto the periodic grid (an atom near the boundary
    //    wraps), and the box ends where the product Gaussian < kScreenEps (a smooth tail -> no ringing).
    // Reuses the PUBLIC Gaussian data (exponents/coeffs/center + pol + norm); the product exponent/center size
    // the box (contracted radials: sized by the most-diffuse primitive pair, value is the full product).

    // Enumerate the screened cross-cell offsets R for the ordered pair (i,j): chi_i (home, at R_i) vs chi_j at
    // R_j+R.  Included when the centres are within reach_i+reach_j (reach = sqrt(-ln eps / alpha_min)); this is
    // the SAME magnitude screen the 1E lattice sums use (consistency = correctness for H psi = eps S psi).
    // The callback receives the INTEGER cell index n (for the caller-supplied Bloch phase e^{ik.R_n} -- the
    // k-convention stays lattice-side) plus its Cartesian translation Roff.
    template <class F>
    void ForImageOffsets(size_t i, size_t j, const UnitCell& A, F&& cb) const
    {
        const rvec3_t Ri=radials[i]->GetCenter(), Rj=radials[j]->GetCenter();
        double aMinI=radials[i]->GetExponents()[0]; for (double e:radials[i]->GetExponents()) aMinI=std::min(aMinI,e);
        double aMinJ=radials[j]->GetExponents()[0]; for (double e:radials[j]->GetExponents()) aMinJ=std::min(aMinJ,e);
        // Screen on the product PREFACTOR exp(-aMinI aMinJ/(aMinI+aMinJ) |Delta|^2) < eps: the pair overlap
        // decays with the centre separation |Delta|=|R_i-(R_j+Roff)|, so include Roff only within this radius.
        const double rr=std::sqrt(-std::log(kScreenEps)*(1.0/aMinI+1.0/aMinJ));
        const rvec3_t dij=Ri-Rj;
        for (const auto& n : A.CellsInSphere(rr+A.GetMaximumCellEdge()))
        {
            const rvec3_t Roff=A.ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z)));
            const rvec3_t d=dij-Roff;                                 // R_i - (R_j + Roff): the centre separation
            if (d.x*d.x+d.y*d.y+d.z*d.z <= rr*rr) cb(n, Roff);
        }
    }
    // The pair's TIGHTEST product resolution: per-component tightest primitive (contracted radials carry a range).
    double MaxExponent(size_t i) const
    {
        double a=radials[i]->GetExponents()[0];
        for (double e:radials[i]->GetExponents()) a=std::max(a,e);
        return a;
    }
    // REL_CUTOFF pair->level assignment (CP2K gaussian_gridlevel): the COARSEST level (largest l) whose cutoff
    // still resolves the pair's TIGHTEST product term alpha_i+alpha_j, req = kRelSafety x the fine grid's own
    // relative resolution ecut_fine (alpha_i+alpha_j)/(2 alpha_max).  kRelSafety=2: the CHARGE-calibrated fine
    // ratio (ecut/p ~ 5 for the auto floor) leaves ~e^{-ecut/2p}=e^{-2.5} spectral tails on a pair at its own
    // cutoff -- fine for the charge (~1e-5 loss) but a few-mHa ENERGY term once V_H/V_xc couple to it (Si SR
    // Gamma sat ~5 mHa below CP2K with kRelSafety=1); doubling the requirement moves only the MID pairs one
    // level up (the diffuse tail stays coarse) and squares the tail to ~e^{-5}.  CP2K's REL_CUTOFF default
    // (30 Ry) corresponds to a ~3x stiffer ratio than the auto floor, so 2 is still on the lean side.
    // For an uncontracted basis this IS the per-primitive-product assignment.
    static constexpr double kRelSafety=2.0;
    size_t PairLevel(size_t i, size_t j, const std::vector<double>& ecut_L, double relCutoffScale) const
    {
        const double req=relCutoffScale*kRelSafety*ecut_L[0]*(MaxExponent(i)+MaxExponent(j))/(2.0*MaxExponent());
        size_t L=0;
        for (size_t l=0; l<ecut_L.size(); l++) if (ecut_L[l]>=req) L=l;   // ecut_L descending (finest first)
        return L;
    }
    // Iterate the (i, j@Roff) product's compact exp-tail box on the N-division grid of cell A, calling
    // f(raster_index, chi_i(r) chi_j(r-R_j-Roff)) at each screened-in, MODULO-WRAPPED grid point.  Shared by
    // collocation (scatter) and integrate-back (gather) -> exact adjoints (same box, chi eval, wrap).
    // THREE screens keep this fast; all are conservative (they drop only sub-eps values):
    //  (1) PREFACTOR-SHRUNK reach: the slowest product term is exp(-p|r-P|^2) x exp(-aI aJ/p |Ri-Rj|^2), so for
    //      a far cross-cell image the exp-tail radius shrinks toward zero -- a screened-in offset near the
    //      reach edge costs a near-empty box instead of the full home-pair box (the dominant win: the
    //      cross-cell offset sum used to re-pay the full box volume per offset);
    //  (2) an exact per-axis FRACTIONAL bounding box of the Cartesian reach-sphere, reach*||row(A^-1)|| (a
    //      sphere spans MORE fractional width than reach/minEdge in a skewed cell -- the old minEdge formula
    //      under-covered an FCC cell by sqrt(3)/sqrt(2), a ~1e-7 clipped charge tail);
    //  (3) an ELLIPSOID pre-screen per point on the slowest-decaying exponent, skipping the exp/poly evals on
    //      the box corners (the exp calls are the kernel's unit of cost; margin e^12 >> any poly growth here).
    //! \a epsEff: the value tolerance this box must honour (default = the geometry screen).  The D-aware
    //! on-the-fly path passes max(kScreenEps, kDensityEps/|weight|) -- a small-weight pair keeps a smaller
    //! box (or none: the prefactor early-out below is the whole-term kill).
    template <class F>
    void ForPairBox(size_t i, size_t j, const rvec3_t& Roff, const UnitCell& A, const ivec3_t& N, F&& f,
                    double epsEff=kScreenEps) const
    {
        const rvec_t ei=radials[i]->GetExponents(), gi=radials[i]->GetCoeffs();
        const rvec_t ej=radials[j]->GetExponents(), gj=radials[j]->GetCoeffs();
        const rvec3_t Ri=radials[i]->GetCenter(), Rj=radials[j]->GetCenter()+Roff;   // chi_j imaged to cell Roff
        const double ni=ns[i], nj=ns[j];
        double aMinI=ei[0]; for (double e:ei) aMinI=std::min(aMinI,e);
        double aMinJ=ej[0]; for (double e:ej) aMinJ=std::min(aMinJ,e);
        const double pMin=aMinI+aMinJ;                       // slowest-decaying product term -> the box size
        const rvec3_t P=(aMinI*Ri+aMinJ*Rj)/pMin;            // its product centre (box centre)
        const double lnE=-std::log(epsEff);
        const rvec3_t dij=Ri-Rj;
        const double pfExp=(aMinI*aMinJ/pMin)*(dij.x*dij.x+dij.y*dij.y+dij.z*dij.z); // prefactor exponent (screen 1)
        if (pfExp>=lnE) return;                              // the whole box is below eps
        const double reach=std::sqrt((lnE-pfExp)/pMin)+1.0;  // shrunk exp-tail radius (+1 a.u. poly margin)
        const rvec3_t fP=A.ToFractional(P);
        // Exact fractional bounding box of the reach-sphere (screen 2): hw_frac(axis) = reach*||(A^-1) row axis||
        // (ToFractional(e_b) is column b of A^-1, so row norms assemble from the three unit-vector images).
        const rvec3_t fex=A.ToFractional(rvec3_t(1,0,0)), fey=A.ToFractional(rvec3_t(0,1,0)), fez=A.ToFractional(rvec3_t(0,0,1));
        const double rbx=std::sqrt(fex.x*fex.x+fey.x*fey.x+fez.x*fez.x);
        const double rby=std::sqrt(fex.y*fex.y+fey.y*fey.y+fez.y*fez.y);
        const double rbz=std::sqrt(fex.z*fex.z+fey.z*fey.z+fez.z*fez.z);
        const int hwx=int(std::ceil(reach*rbx*N.x))+1, hwy=int(std::ceil(reach*rby*N.y))+1, hwz=int(std::ceil(reach*rbz*N.z))+1;
        const long cx=std::lround(fP.x*N.x), cy=std::lround(fP.y*N.y), cz=std::lround(fP.z*N.z);
        const double lnCut=lnE+12.0;                         // ellipsoid pre-screen bound (screen 3)
        // INCREMENTAL grid walk: r = A(g/N) advances by a CONSTANT lattice step per index, so the per-point
        // ToCartesian call (a 3x3 matrix multiply + PLT hop; ~14% of the box-loop profile) hoists to three
        // axis steps + adds.  Accumulation drift over <=few-hundred steps is ~1e-13 -- far below kScreenEps.
        const rvec3_t sx=A.ToCartesian(rvec3_t(1.0/N.x,0,0)), sy=A.ToCartesian(rvec3_t(0,1.0/N.y,0)),
                      sz=A.ToCartesian(rvec3_t(0,0,1.0/N.z));
        const rvec3_t r00=A.ToCartesian(rvec3_t(double(cx-hwx)/N.x, double(cy-hwy)/N.y, double(cz-hwz)/N.z));
        rvec3_t rx=r00;
        for (int dx=-hwx; dx<=hwx; dx++, rx=rx+sx)
        {
            rvec3_t ry=rx;
            for (int dy=-hwy; dy<=hwy; dy++, ry=ry+sy)
            {
                rvec3_t r=ry;
                for (int dz=-hwz; dz<=hwz; dz++, r=r+sz)
                {
                    const rvec3_t di=r-Ri, dj=r-Rj;
                    const double ri2=di.x*di.x+di.y*di.y+di.z*di.z, rj2=dj.x*dj.x+dj.y*dj.y+dj.z*dj.z;
                    if (aMinI*ri2+aMinJ*rj2 > lnCut) continue;   // slowest term < eps*e^-12 -> skip exp/poly evals
                    double radI=0.0; for (size_t p=0;p<ei.size();p++) radI+=gi[p]*std::exp(-ei[p]*ri2);
                    double radJ=0.0; for (size_t q=0;q<ej.size();q++) radJ+=gj[q]*std::exp(-ej[q]*rj2);
                    const double val=(ni*pols[i](di)*radI)*(nj*pols[j](dj)*radJ);
                    if (std::fabs(val)<epsEff) continue;
                    const long gx=cx+dx, gy=cy+dy, gz=cz+dz;     // unwrapped grid index (may lie outside [0,N))
                    const long mx=((gx%N.x)+N.x)%N.x, my=((gy%N.y)+N.y)%N.y, mz=((gz%N.z)+N.z)%N.z;  // wrap
                    f((size_t(mx)*N.y+my)*N.z+mz, val);
                }
            }
        }
    }
    // --- STREAM CACHE for the analytic collocation / integrate-back ----------------------------------------
    // The per-(pair, cross-cell offset) box streams -- wrapped grid index + analytic pair value -- are pure
    // GEOMETRY: identical across SCF iterations AND across k-blocks (the Bloch phase enters only at the
    // contraction), while their evaluation (exp/poly per point) is essentially 100% of the GPW DFT tier's
    // per-iteration cost (perf 2026-07-14: exp 24%, box-loop bodies 41%, ToCartesian 14%...).  So each stream
    // is built ONCE per grid-ladder shape and replayed as a pure gather/scatter -- the analytic sibling of the
    // deleted PhiOnGrid cache, but per-pair-compact instead of dense.  One cache per distinct ladder shape
    // (the SCF ladder and the K=1 local-PP / unit-gate shape coexist); replay order == build order == the
    // uncached loop order, so results are BIT-IDENTICAL to the uncached path.
    //! TWO value tiers: fp64 (\c val) for pairs inside the primary budget, fp32 (\c val32) for the overflow
    //! tier -- exactly one of the two is filled per stream.  \c maxv (the stream's largest |value|) powers the
    //! D-AWARE replay kill: an offset whose |weight|*maxv is below the density tolerance contributes nothing
    //! resolvable and is skipped whole.
    struct PairOffsetStream { ivec3_t n; std::vector<unsigned> idx; std::vector<double> val; std::vector<float> val32; double maxv=0.0; };
    struct PairStreams     { size_t level=0; bool cached=false; std::vector<PairOffsetStream> offsets; };
    struct StreamCache
    {
        std::vector<ivec3_t> N_L; std::vector<double> ecut_L;   // the ladder shape (snapshot key)
        double relCutoffScale=1.0;                              //   + the field-sharpness assignment scale
        std::vector<PairStreams> pairs;                         // indexed [i*n+j], j>=i only
    };
    static constexpr size_t kMaxStreamCaches=4;                 // ladder + K=1 + unit-gate shapes; guards churn
    //! Global stream-cache budgets in POINTS, TWO TIERS.  A diffuse basis in a small cell keeps hundreds of
    //! screened offsets per pair -- caching them ALL blew past physical RAM on NaF (uncapped: 27 GB virt on a
    //! 16 GB box; measured demand 952M pts).  Tier 1 (fp64, 12 B/pt, ~1.8 GB): replay is BIT-IDENTICAL to
    //! on-the-fly evaluation; sized so the Si SR ladder caches completely (demand 104.9M) -- every Si anchor
    //! and machine-precision kernel gate lives entirely in this tier.  Tier 2 (fp32 values, 8 B/pt, ~5.6 GB):
    //! the overflow pairs store float values instead of falling to on-the-fly -- replay noise is ~6e-8
    //! RELATIVE per value (NaF anchors are 1e-3-scale; Tr(DS) charge is analytic, unaffected), and the
    //! collocate/integrate ADJOINT stays machine-exact because both directions replay the SAME stream.
    //! Pairs beyond BOTH budgets fall back to on-the-fly evaluation (correct, slower).
    static constexpr size_t kStreamBudgetPts   =150'000'000;   // tier 1: fp64, bit-identical replay
    //! Tier 2 sized so NaF's MEASURED demand (952M total: 150M fp64 + 802M fp32) caches COMPLETELY -- at
    //! 700M its last 314 (small) pairs re-evaluated on the fly every sweep.  Peak stream RAM ~ 8.6 GB.
    static constexpr size_t kStreamBudgetPtsF32=850'000'000;   // tier 2: fp32 values (the NaF coverage lever)
    const StreamCache& EnsureStreams(const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                     const std::vector<double>& ecut_L, double relCutoffScale) const
    {
        for (const auto& c : itsStreamCaches)
        {
            if (c.relCutoffScale!=relCutoffScale || c.ecut_L!=ecut_L || c.N_L.size()!=N_L.size()) continue;
            bool same=true;
            for (size_t l=0;l<N_L.size();l++)
                if (c.N_L[l].x!=N_L[l].x || c.N_L[l].y!=N_L[l].y || c.N_L[l].z!=N_L[l].z) { same=false; break; }
            if (same) return c;
        }
        if (itsStreamCaches.size()>=kMaxStreamCaches) itsStreamCaches.clear();   // unexpected shape churn
        itsStreamCaches.emplace_back();
        StreamCache& c=itsStreamCaches.back();
        c.N_L=N_L; c.ecut_L=ecut_L; c.relCutoffScale=relCutoffScale;
        const size_t n=size();
        c.pairs.resize(n*n);
        size_t budget64=kStreamBudgetPts, budget32=kStreamBudgetPtsF32;
        for (size_t cc=0; cc<itsStreamCaches.size()-1; cc++)          // budgets are GLOBAL across cache shapes
            for (const auto& ps : itsStreamCaches[cc].pairs)
                for (const auto& st : ps.offsets)
                    if (!st.val.empty()) budget64 -= std::min(budget64, st.idx.size());
                    else                 budget32 -= std::min(budget32, st.idx.size());
        size_t nPairs=0, nCached64=0, nCached32=0, pts64=0, pts32=0, ptsDropped=0;
        for (auto i:indices()) for (auto j:indices(i))
        {
            PairStreams& ps=c.pairs[i*n+j];
            ps.level=PairLevel(i,j,ecut_L,relCutoffScale);
            nPairs++;
            if (budget64==0 && budget32==0) { ptsDropped++; continue; }   // both tiers exhausted
            const ivec3_t N=N_L[ps.level];
            size_t pts=0;
            ForImageOffsets(i,j,A,[&](const ivec3_t& nn, const rvec3_t& Roff)
            {
                PairOffsetStream st; st.n=nn;
                ForPairBox(i,j,Roff,A,N,[&](size_t idx,double v)
                           { st.idx.push_back(unsigned(idx)); st.val.push_back(v);
                             st.maxv=std::max(st.maxv,std::fabs(v)); });
                pts+=st.idx.size();
                if (!st.idx.empty()) ps.offsets.push_back(std::move(st));
            });
            // Tier the pair: fp64 while the primary budget holds (bit-identical replay), fp32 for the
            // overflow tier, drop only past BOTH.  Per-pair skip-and-continue, never a lockout -- a later
            // (smaller) pair may still fit a tier.  (The old budget=0 lockout un-cached every subsequent
            // pair after the first oversized one -- the dominant cost in the 2026-07-15 multi-k profile.)
            if (pts<=budget64)
            {
                budget64-=pts; pts64+=pts; nCached64++;
                ps.cached=true;
            }
            else if (pts<=budget32)
            {
                for (auto& st : ps.offsets)                       // demote values to fp32 (idx stays exact)
                {
                    st.val32.assign(st.val.begin(), st.val.end());
                    st.val.clear(); st.val.shrink_to_fit();
                }
                budget32-=pts; pts32+=pts; nCached32++;
                ps.cached=true;
            }
            else { ps.offsets.clear(); ps.offsets.shrink_to_fit(); ptsDropped+=pts; }
        }
        // One-line readout per cache build (static setup, not per-iteration): the budget headroom is THE lever
        // for the analytic path's speed, so make coverage visible (dropped pairs re-evaluate every iteration).
        std::cerr << "[stream cache] shape=(";
        for (size_t l=0;l<N_L.size();l++) std::cerr << (l?",":"") << N_L[l].x;
        std::cerr << ") scale=" << relCutoffScale << "  pairs " << nPairs
                  << ": fp64 " << nCached64 << " (" << pts64 << " pts), fp32 " << nCached32
                  << " (" << pts32 << " pts), dropped " << ptsDropped
                  << " pts (budgets " << kStreamBudgetPts << "/" << kStreamBudgetPtsF32 << ")" << std::endl;
        return c;
    }

    // Collocate the grid density onto the multi-grid ladder (Molecule::LatticeSum1E::CollocateDensity).  The
    // cell-periodic density of the Bloch orbital products chi_i^k conj(chi_j^k) is, per (pair, cross-cell
    // offset R):  rho(r) = Sum_ij Sum_R Re[D_ij e^{-ik.R}] chi_i(r) chi_j(r-R)   (chi real Gaussians; the
    // modulo wrap carries the cell tiling).  The offset phase is the CONJUGATE e^{-ik.R} -- the density
    // convention (rho = Sum D_ij chi_i^k conj(chi_j^k); the ket slot conjugated, doc/GPWPlan.md complex-k fix);
    // summing the full n^2 with the Re[] weight is exact (the (j,i,-R) partner term supplies the conjugate).
    // Each pair scatters on ITS level (PairLevel); Integral of the total = Tr(D S^k) to screening tolerance.
    std::vector<rvec_t> CollocateDensity(const chmat_t& D, const cellphase_t& phase, const UnitCell& A,
                                         const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L) const
    {
        const size_t K=N_L.size();
        assert(K>0 && ecut_L.size()==K);
        std::vector<rvec_t> rho(K);
        for (size_t l=0; l<K; l++) rho[l]=rvec_t(size_t(N_L[l].x)*N_L[l].y*N_L[l].z, 0.0);
        // Hermitian fold: the (j,i,-R) term contributes the SAME (idx, value, weight) as (i,j,R) -- the product
        // field is the lattice-translated twin and D_ji e^{-ik(-R)} = conj(D_ij e^{-ikR}) -- so loop j>=i and
        // double the off-diagonal weight (a 2x saving over the full n^2 loop, exact).
        const StreamCache& sc=EnsureStreams(A,N_L,ecut_L,1.0);   // density: the smooth-field calibration
        const size_t nn=size();
        for (auto i:indices()) for (auto j:indices(i))
        {
            const dcmplx Dij(D(i,j));
            if (Dij==dcmplx(0.0)) continue;
            const double fold=(i==j)?1.0:2.0;
            const PairStreams& ps=sc.pairs[i*nn+j];
            rvec_t& r=rho[ps.level];
            if (ps.cached)
                for (const PairOffsetStream& st : ps.offsets)
                {
                    const double c=fold*std::real(Dij*std::conj(phase(st.n)));  // Re[D_ij e^{-ik.R_n}] offset weight
                    if (std::fabs(c)*st.maxv < kDensityEps) continue;   // D-aware kill: sub-eps density term
                    const unsigned* ix=st.idx.data();
                    if (!st.val.empty())
                    {
                        const double* v=st.val.data();
                        for (size_t k=0, m=st.idx.size(); k<m; k++) r[ix[k]]+=c*v[k];
                    }
                    else                                            // fp32 overflow tier (values demoted)
                    {
                        const float* v=st.val32.data();
                        for (size_t k=0, m=st.idx.size(); k<m; k++) r[ix[k]]+=c*double(v[k]);
                    }
                }
            else                                                    // over the cache budget: evaluate on the fly
                ForImageOffsets(i,j,A,[&](const ivec3_t& n, const rvec3_t& Roff)
                {
                    const double c=fold*std::real(Dij*std::conj(phase(n)));
                    if (c==0.0) return;
                    // D-aware continuous shrink: this box only needs accuracy kDensityEps/|c| (clamped -- a
                    // |c|>1 never grows past the geometry screen); the prefactor early-out inside is the kill.
                    ForPairBox(i,j,Roff,A,N_L[ps.level],[&](size_t idx,double v){r[idx]+=c*v;},
                               std::max(kScreenEps, kDensityEps/std::fabs(c)));
                });
        }
        return rho;
    }
    // --- Phase-independent integrate-back memo ---------------------------------------------------------------
    // h_ij(k) = w_l Sum_n e^{+ik.R_n} B_ij(n), with B_ij(n) = Sum_box chi_i chi_j^n V.  The expensive per-offset
    // reductions B are k-INDEPENDENT (the Bloch phase enters only the final contraction), and the SAME field V
    // is integrated repeatedly: the static local PP by EVERY k-block (once per SCF each), and the per-iteration
    // KS fields once per k-block within an iteration.  So memoize {B_ij(n)} keyed on the EXACT (ladder shape,
    // scale, V_L): the first caller pays the sweep, every other k-block only contracts phases (2026-07-15
    // multi-k profile: the per-k MakeLocalPP sweep alone was ~10% of the anchor).  Like itsStreamCaches, the
    // cell A is assumed fixed per basis instance (the centres ARE this basis's data) and is not in the key.
    // Field equality is EXACT per-element (never blaze relaxed equal); contraction order == the direct
    // evaluation's offset order, so a memo hit is BIT-IDENTICAL to recomputation.
    struct PairB         { size_t level=0; std::vector<std::pair<ivec3_t,double>> nb; };  // (offset, B_ij(n))
    struct IntegrateMemo
    {
        std::vector<ivec3_t> N_L; std::vector<double> ecut_L;   // the ladder shape (snapshot key)
        double relCutoffScale=1.0;                              //   + the field-sharpness assignment scale
        std::vector<rvec_t>  V_L;                               //   + the integrated field itself (exact)
        std::vector<PairB>   B;                                 // indexed [i*n+j], j>=i only
    };
    static constexpr size_t kMaxIntegrateMemos=4;               // static PP + the iteration's KS fields
    bool SameShape(const IntegrateMemo& m, const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L,
                   double relCutoffScale) const
    {
        if (m.relCutoffScale!=relCutoffScale || m.ecut_L!=ecut_L || m.N_L.size()!=N_L.size()) return false;
        for (size_t l=0;l<N_L.size();l++)
            if (m.N_L[l].x!=N_L[l].x || m.N_L[l].y!=N_L[l].y || m.N_L[l].z!=N_L[l].z) return false;
        return true;
    }
    static bool SameField(const std::vector<rvec_t>& a, const std::vector<rvec_t>& b)
    {
        if (a.size()!=b.size()) return false;
        for (size_t l=0;l<a.size();l++)
        {
            if (a[l].size()!=b[l].size()) return false;
            for (size_t k=0, m=a[l].size(); k<m; k++) if (a[l][k]!=b[l][k]) return false;   // exact, not blaze equal
        }
        return true;
    }
    // Integrate-back (the collocation ADJOINT): h_ij = <chi_i^k|V|chi_j^k> = Sum_R e^{+ik.R} w_l Sum_box
    // chi_i chi_j^R V, over the SAME screened offsets + compact boxes + wrap + level assignment as collocation
    // -> exact adjoint (Integral rho.V == Tr(D h) to machine precision, variational).  The offset phase here is
    // the DIRECT e^{+ik.R} (the Bloch law on the ket image; conjugate of the density side's weight -- the same
    // +/- pairing as the KB projector's conj(phase), doc/GPWPlan.md complex-k fix).  Hermitian; real at Gamma.
    // Only V is sampled (weighted by the analytic Gaussians), never the sharp orbital product.
    //! \a screenD: OPTIONAL density-magnitude screen (the D-AWARE ACTIVE SET).  When supplied, each
    //! (pair, offset) is kept exactly when the DENSITY-side collocation of \a screenD keeps it
    //! (|weight|*max|value| >= kDensityEps, the identical criterion) -- so within an SCF iteration the
    //! collocate/integrate ADJOINT is machine-exact on the shared truncated operator, and the sweep only
    //! touches terms the density can resolve.  Screened calls bypass the memo (they are cheap by
    //! construction, and the active set changes with D while the memo is keyed on V alone).
    chmat_t IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                               const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L,
                               double relCutoffScale=1.0, const chmat_t* screenD=nullptr) const
    {
        const size_t K=N_L.size();
        assert(K>0 && ecut_L.size()==K && V_L.size()==K);
        const size_t nn=size();
        chmat_t h(size());
        const bool memoize = (screenD==nullptr);
        // Memo hit: the per-offset reductions are already computed for this exact field -- contract phases only.
        if (memoize)
            for (const auto& m : itsIntegrateMemos)
                if (SameShape(m,N_L,ecut_L,relCutoffScale) && SameField(m.V_L,V_L))
                {
                    for (auto i:indices()) for (auto j:indices(i))
                    {
                        const PairB& pb=m.B[i*nn+j];
                        const double w=A.GetCellVolume()/double(V_L[pb.level].size());
                        dcmplx s(0.0);
                        for (const auto& [n,b] : pb.nb) s+=phase(n)*b;
                        s*=w;
                        h(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;
                    }
                    return h;
                }
        // Miss: compute (streams replay where cached, on-the-fly otherwise), recording B for the next caller.
        // relCutoffScale != 1 marks a STATIC sharp-field call (the local PP, built once per SCF): evaluate on
        // the fly -- caching a single-shot sweep wastes the whole budget the per-iteration path needs.  (Its
        // B reductions land in the memo, so the other k-blocks never repeat the sweep.)
        const StreamCache* sc = (relCutoffScale==1.0) ? &EnsureStreams(A,N_L,ecut_L,relCutoffScale) : nullptr;
        IntegrateMemo* memo=nullptr;
        if (memoize)
        {
            if (itsIntegrateMemos.size()>=kMaxIntegrateMemos) itsIntegrateMemos.erase(itsIntegrateMemos.begin());
            itsIntegrateMemos.emplace_back();
            memo=&itsIntegrateMemos.back();
            memo->N_L=N_L; memo->ecut_L=ecut_L; memo->relCutoffScale=relCutoffScale; memo->V_L=V_L;
            memo->B.resize(nn*nn);
        }
        for (auto i:indices()) for (auto j:indices(i))       // j>=i (Hermitian upper triangle)
        {
            const size_t  l = sc ? sc->pairs[i*nn+j].level : PairLevel(i,j,ecut_L,relCutoffScale);
            const rvec_t& V=V_L[l];
            const double  w=A.GetCellVolume()/double(V.size());   // the level's quadrature weight Omega/Npts(l)
            PairB  dummy;
            PairB& pb = memo ? memo->B[i*nn+j] : dummy;
            pb.level=l;
            // The D-aware offset weight: the SAME |c| the density collocation of screenD would apply --
            // fold*Re[D_ij e^{-ik.R_n}] -- so both directions keep the identical active set.
            const double fold=(i==j)?1.0:2.0;
            const dcmplx Dij = screenD ? dcmplx((*screenD)(i,j)) : dcmplx(0.0);
            dcmplx s(0.0);
            if (sc && sc->pairs[i*nn+j].cached)
                for (const PairOffsetStream& st : sc->pairs[i*nn+j].offsets)
                {
                    if (screenD &&
                        std::fabs(fold*std::real(Dij*std::conj(phase(st.n))))*st.maxv < kDensityEps) continue;
                    double b=0.0;
                    const unsigned* ix=st.idx.data();
                    if (!st.val.empty())
                    {
                        const double* v=st.val.data();
                        for (size_t k=0, m=st.idx.size(); k<m; k++) b+=v[k]*V[ix[k]];
                    }
                    else                                            // fp32 overflow tier (values demoted)
                    {
                        const float* v=st.val32.data();
                        for (size_t k=0, m=st.idx.size(); k<m; k++) b+=double(v[k])*V[ix[k]];
                    }
                    pb.nb.emplace_back(st.n,b);
                    s+=phase(st.n)*b;
                }
            else                                                    // static sharp-field call or over-budget pair
                ForImageOffsets(i,j,A,[&](const ivec3_t& n, const rvec3_t& Roff)
                {
                    double epsEff=kScreenEps;
                    if (screenD)
                    {
                        const double c=std::fabs(fold*std::real(Dij*std::conj(phase(n))));
                        if (c==0.0) return;
                        epsEff=std::max(kScreenEps, kDensityEps/c);
                    }
                    double b=0.0;
                    ForPairBox(i,j,Roff,A,N_L[l],[&](size_t idx,double v){b+=v*V[idx];}, epsEff);
                    pb.nb.emplace_back(n,b);
                    s+=phase(n)*b;
                });
            s*=w;
            h(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) auto-set to conj
        }
        return h;
    }

    // --- 3-centre (DFT) and 4-centre (HF) kernels ---------------------------------------------------
    // General multi-evaluator elements: each (evaluator, index) pair names one basis component, so the
    // same kernel serves both Coulomb/Exchange (the caller just maps its loop indices into the slots).
    // `this` is the A slot.  Normalizations of every slot are folded in, matching MakeOverlap3C /
    // MakeDirect / MakeExchange.  (A member may read another same-type evaluator's private data.)

    // <ab|c> 3-centre integrals (M&D, evaluated on the C radial) -- one named function per kind.
    double OverlapThreeC(size_t iA, const NR_Evaluator& eB, size_t iB, const NR_Evaluator& eC, size_t iC) const
    {
        return eC.radials[iC]->Overlap3C(*radials[iA], *eB.radials[iB], pols[iA], eB.pols[iB], eC.pols[iC])
             * ns[iA] * eB.ns[iB] * eC.ns[iC];
    }
    double RepulsionThreeC(size_t iA, const NR_Evaluator& eB, size_t iB, const NR_Evaluator& eC, size_t iC) const
    {
        return eC.radials[iC]->Repulsion3C(*radials[iA], *eB.radials[iB], pols[iA], eB.pols[iB], eC.pols[iC])
             * ns[iA] * eB.ns[iB] * eC.ns[iC];
    }

    // (ab|cd) 4-centre electron-repulsion integral (M&D, evaluated on the D radial).
    double FourC(size_t iA,
                 const NR_Evaluator& eB, size_t iB,
                 const NR_Evaluator& eC, size_t iC,
                 const NR_Evaluator& eD, size_t iD) const
    {
        return eD.radials[iD]->Repulsion4C(*radials[iA], *eB.radials[iB], *eC.radials[iC],
                                           pols[iA], eB.pols[iB], eC.pols[iC], eD.pols[iD])
             * ns[iA] * eB.ns[iB] * eC.ns[iC] * eD.ns[iD];
    }


private:
    mutable std::vector<StreamCache>   itsStreamCaches;    //!< analytic pair-box streams, one per ladder shape
    mutable std::vector<IntegrateMemo> itsIntegrateMemos;  //!< phase-independent B_ij(n) per exact field (LRU)
};

static_assert(is1E_Evaluator <NR_Evaluator>, "NR_Evaluator must satisfy is1E_Evaluator");
static_assert(isDFT_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isDFT_Evaluator (ThreeC kernel)");
static_assert( isHF_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isHF_Evaluator (FourC kernel)");

} //namespace

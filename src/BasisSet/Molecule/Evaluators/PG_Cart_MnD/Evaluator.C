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
#include <fstream>   // /proc/self/statm (GPW_RSS_TRACE diagnostics)
#include <iostream>   // std::cerr (the one-line stream-cache coverage readout)
#include <vector>
#include <cmath>      // std::sqrt/std::log (the lattice-sum magnitude-screening reach radius)
#include <algorithm>  // std::min (alpha_min per radial)
#include <functional> // cellphase_t (the caller-supplied Bloch phase of a cross-cell offset)
#include <utility>    // std::pair (the flattened (i,j) pair list for the OpenMP loop)
#include <map>        // the stream fold's offset -> orbit-multiplicity tables (T3 route (b))
#include <optional>   // the conditionally-charged per-iteration timing bucket (integrate-back)
#include <exception>  // std::exception_ptr (throw containment across the OpenMP pair loops)
#include <memory>     // std::shared_ptr (the per-atom operator GaussianRF, shared across its polynomial terms)
#include <cstdlib>    // std::getenv/std::atoi (the GPW_OMP_THREADS opt-in knob)
// GPW pair-loop parallelism (opt-in via GPW_OMP_THREADS).  Gated on QCHEM_OPENMP -- our OWN macro (defined
// by CMake when the OpenMP flags are applied) rather than _OPENMP, because this LLVM toolchain ships no
// libomp and the libgomp fallback (-fopenmp=libgomp) honours the pragmas but does NOT define _OPENMP.  The
// private-buffer + critical-reduce pattern below needs no <omp.h> (no omp_*() calls), only the pragmas.
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;      // PGData
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;  // GaussianRF named kernels
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;// Polarization
import qchem.BasisSet.Molecule.LatticeSum1E;                       // GaussianFunction (the <chi_i|g> seam type)
export import qchem.Symmetry.Lattice_3D.SpaceGroup;                // DirectOp {W|τ} (the T3 stream-fold ops)
import qchem.Structure;
import qchem.UnitCell;                                             // UnitCell (ToCartesian/ToFractional: grid<->cell for collocation)
import qchem.Types;
import qchem.Blaze;                                                // rsmat_t (the lattice-sum matrices)
import qchem.Reporting;                                            // fold the stream-cache coverage into grids.stream

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
    // Magnitude screen (CP2K EPS_PGF_ORB analog).  Default 1e-10 (numerically exact for GPW).
    // Env instrument GPW_SCREEN_EPS raises it to SHRINK the lattice-sum reach + collocation boxes
    // for diffuse bases (demo/robustness runs -- drops only sub-eps terms; not for production digits).
    static double kScreenEps() { static const double e=[]{const char* s=std::getenv("GPW_SCREEN_EPS"); return s?std::atof(s):1e-10;}(); return e; }
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
    //! D-aware density-magnitude screen (sizes the collocation boxes).  Env instrument GPW_DENSITY_EPS.
    static double kDensityEps() { static const double e=[]{const char* s=std::getenv("GPW_DENSITY_EPS"); return s?std::atof(s):1e-10;}(); return e; }

    //! The lattice-sum ECONOMY report (LatticeSum1E face; doc there).  Reach numbers use the SAME formulas
    //! as the screens: a single orbital's tail reach is sqrt(-ln eps/alpha), a pair's conservative reach is
    //! the sum -- worst pair = diffuse x diffuse = 2*sqrt(-ln eps/alpha_min).  The CellsInSphere count over
    //! that reach is the per-pair image-enumeration size a grad student sees JUMP when diffuse functions
    //! arrive; the screens then prune it per term (kept counts = the [stream cache] offsets readout).
    void EmitLatticeSumReport(const UnitCell& A) const
    {
        const double aMin=MinExponent(), aMax=MaxExponent();
        const double reachS =2.0*std::sqrt(-std::log(kScreenEps ())/aMin);   // analytic 1E/V_local enumeration
        const double reachD =2.0*std::sqrt(-std::log(kDensityEps())/aMin);   // collocation offset enumeration
        const long   cellsS =(long)A.CellsInSphere(reachS).size();
        const long   cellsD =(long)A.CellsInSphere(reachD).size();
        std::cout<<"[lattice sums] alpha_min="<<aMin<<" alpha_max="<<aMax
                 <<"  analytic 1E/V_local: eps="<<kScreenEps()<<" (GPW_SCREEN_EPS) pair reach="<<reachS
                 <<" au = "<<cellsS<<" cells;  collocation offsets: eps="<<kDensityEps()
                 <<" (GPW_DENSITY_EPS) pair reach="<<reachD<<" au = "<<cellsD
                 <<" cells (kept counts on the [stream cache] line)"<<std::endl;
        qchem::report::EmitAt("grids", "latticeSums", {
            {"alphaMin",aMin}, {"alphaMax",aMax},
            {"screenEps",kScreenEps()},  {"pairReach1E",reachS},   {"cells1E",cellsS},
            {"densityEps",kDensityEps()},{"pairReachColl",reachD}, {"cellsColl",cellsD}});
    }
    //! The Bloch phase of an integer cell offset (== Molecule::LatticeSum1E::cellphase_t -- the same
    //! std::function type; the k-CONVENTION stays with the lattice-side caller, Gamma = the constant 1).
    using cellphase_t = std::function<dcmplx(const ivec3_t& n)>;
    // The series is summed to eps INTERNALLY per shell pair (ForImageOffsets -- the SAME exact-threshold
    // magnitude screen the collocation kernels use): THERE IS NO CUT in the R direction (doc/GPWPlan.md
    // pin).  The caller supplies only the phase oracle + the cell; no radius or translation list exists.
    // Hermitian: fill the upper triangle (the (i,j) Bloch element); the (i,i) enumeration is inversion-
    // symmetric so the diagonal is real (projected explicitly against roundoff).
    template <class Kernel> chmat_t LatticeSum(const cellphase_t& phase, const UnitCell& A, Kernel K) const
    {
        chmat_t S(size());
        const bool trace=(std::getenv("GPW_RSS_TRACE")!=nullptr);
        // GEOMETRY-CACHE BOUND (the full-SR OOM, doc/GPWPlan.md): the MnD Omega/RNLM/H3 caches key on
        // per-instance IDs, and AtCenter mints fresh clones per (pair, offset) -- a diffuse basis' ~2k
        // images/pair would grow them to GBs inside ONE matrix build.  The size-0 budget evicts each
        // clone as the next is inserted (O(1) resident); safe here because the 1E/3C kernels hold at
        // most one live borrow per cache.  Replaces the old per-pair ClearGeometryCaches().
        const GeometryCacheBudget geoBudget;
        for (auto i:indices())
        {
            for (auto j:indices(i))
            {
                const rvec3_t cj=radials[j]->GetCenter();
                dcmplx s(0.0);
                ForImageOffsets(i,j,A,[&](const ivec3_t& n, const rvec3_t& Roff)
                                { s += phase(n) * K(i, j, radials[j]->AtCenter(cj+Roff)); });
                s *= ns[i]*ns[j];
                S(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) auto-set to conj
            }
            if (trace)
            {
                std::ifstream f("/proc/self/statm"); size_t vm=0, rs=0; f>>vm>>rs;
                std::cerr<<"[rss row] i="<<i<<" rss="<<(rs*4096/1048576)<<" MB"<<std::endl;
            }
        }
        return S;
    }
    chmat_t MakeOverlap(const cellphase_t& phase, const UnitCell& A) const
    {   return LatticeSum(phase,A,[this](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Overlap2C(cj,pols[i],pols[j]);}); }
    chmat_t MakeKinetic(const cellphase_t& phase, const UnitCell& A) const
    {   return LatticeSum(phase,A,[this](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Grad2   (cj,pols[i],pols[j]);}); }
    chmat_t MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const
    {   return LatticeSum(phase,A,[this,cl](size_t i,size_t j,const GaussianRF& cj){return radials[i]->Nuclear(cj,pols[i],pols[j],cl);}); }

    // Molecule::LatticeSum1E: the analytic 3-centre (Overlap3C) MATRIX of a per-atom LOCAL Gaussian operator --
    // <chi_i^k | Sum_a g_a | chi_j^k>, the local pseudopotential's short-range assembly (doc/GPWPlan.md 0e-PP).
    // Each atom's operator g_a = opForZ(Z_a) sits in the OPERATOR (C) slot of Overlap3C, between chi_i and the
    // lattice-imaged chi_j; a compact operator's image atoms are negligible so only home-cell atoms are placed.
    // The LatticeSum's chi_i-chi_j magnitude screen suffices: a screened-in pair whose product is disjoint from
    // g_a gives Overlap3C ~ 0 (the integrand needs all three overlapping), and a screened-OUT pair (chi_i,chi_j
    // disjoint) has product ~ 0 everywhere, so no significant 3-centre term is dropped.
    chmat_t MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                              const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ) const
    {
        const std::vector<OpTerm> ops=BuildOpTerms(cl,opForZ,&A);   // periodic: operator over its screened images
        // 3-CENTRE MAGNITUDE SCREEN (the localized-PP anchor).  The chi_i-chi_j pair screen (ForImageOffsets)
        // enumerates every image a DIFFUSE pair reaches -- but the local PP g is LOCALIZED at its atom, so a
        // (chi_i, chi_j image) product disjoint from g contributes ~0.  Screen each op by the EXACT 3-Gaussian
        // leading magnitude exp(-(ab|AB|^2 + ac|AC|^2 + bc|BC|^2)/(a+b+c)) (M&D's Eabc), using each radial's
        // MOST DIFFUSE primitive (a=MinExponent -> largest reach -> conservative: drops only sub-eps terms), and
        // SKIP the expensive Overlap3C when it is below eps.  This is the 3-centre analogue of the pair reach
        // sqrt(-ln eps (1/a_i+1/a_j)): a diffuse orbital no longer forces a huge operator sum, because g anchors
        // the product to its atom.  (3-centre kernel under LatticeSum's size-0 GeometryCacheBudget: each surviving
        // triple's Hermite3 block is built, used, and evicted -- O(1) resident.)
        const double lne=-std::log(kScreenEps());
        return LatticeSum(phase,A,[this,&ops,lne](size_t i,size_t j,const GaussianRF& cj)
        {
            const rvec3_t Ri=radials[i]->GetCenter(), Rj=cj.GetCenter();   // chi_i home, chi_j imaged
            const double ai=MinExponent(i), aj=MinExponent(j);            // AtCenter preserves exponents
            const rvec3_t AB=Ri-Rj; const double AB2=AB*AB;
            double s=0.0;
            for (const auto& ot : ops)
            {
                const rvec3_t Rc=ot.gr->GetCenter();
                const double c=ot.gr->GetExponents()[0];                  // op is uncontracted (single g.alpha)
                const rvec3_t AC=Ri-Rc, BC=Rj-Rc;
                if ((ai*aj*AB2 + ai*c*(AC*AC) + aj*c*(BC*BC))/(ai+aj+c) > lne) continue;  // prefactor < eps: skip
                s += ot.c * ot.gr->Overlap3C(*radials[i], cj, pols[i], pols[j], ot.pol);
            }
            return s;
        });
    }
    // The FINITE (home-term-only) <chi_i | Sum_a g_a | chi_j>: the molecule/box limit (no lattice), real symmetric.
    chmat_t MakeLocalGaussian(const Structure* cl,
                              const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ) const
    {
        const std::vector<OpTerm> ops=BuildOpTerms(cl,opForZ,nullptr);   // finite: home atoms only
        const double lne=-std::log(kScreenEps());   // SAME 3-centre screen as the periodic form (finite==lattice)
        chmat_t S(size());
        for (auto i:indices()) for (auto j:indices(i))
        {
            const rvec3_t Ri=radials[i]->GetCenter(), Rj=radials[j]->GetCenter();
            const double ai=MinExponent(i), aj=MinExponent(j);
            const rvec3_t AB=Ri-Rj; const double AB2=AB*AB;
            double s=0.0;
            for (const auto& ot : ops)
            {
                const rvec3_t Rc=ot.gr->GetCenter();
                const double c=ot.gr->GetExponents()[0];
                const rvec3_t AC=Ri-Rc, BC=Rj-Rc;
                if ((ai*aj*AB2 + ai*c*(AC*AC) + aj*c*(BC*BC))/(ai+aj+c) > lne) continue;  // prefactor < eps: skip
                s += ot.c * ot.gr->Overlap3C(*radials[i], *radials[j], pols[i], pols[j], ot.pol);
            }
            S(i,j)=dcmplx(s*ns[i]*ns[j], 0.0);   // finite: real symmetric (upper triangle; chmat_t mirrors)
        }
        return S;
    }

    // Molecule::LatticeSum1E: the lattice-summed overlap of every basis function with ONE caller-supplied
    // Cartesian-Gaussian function g -- b_i = Sum_n phase(n) <chi_i | g(.-C-R_n)>.  g stands in the chi_j
    // slot of the 2-centre kernel: each monomial term is a Polarization on an uncontracted GaussianRF, so
    // the sum reuses the analytic M&D Overlap2C verbatim.  The offsets are enumerated INTERNALLY per
    // (chi_i, g) at the exact pair threshold sqrt(-ln eps (1/alpha_min_i + 1/alpha_g)) -- the same
    // eps-converged-series contract as the pair sums; chi_i carries ns[i], g is integrated RAW.
    cvec_t MakeOverlap(const cellphase_t& phase, const UnitCell& A,
                       const Molecule::LatticeSum1E::GaussianFunction& g) const
    {
        const double lne=-std::log(kScreenEps());
        const int Lg=GaussDegree(g);
        cvec_t b(size(), dcmplx(0.0));
        const GeometryCacheBudget geoBudget;   // size-0: g's image clones stay O(1) (see LatticeSum)
        for (auto i:indices())
        {
            double amin=radials[i]->GetExponents()[0];
            for (double e:radials[i]->GetExponents()) amin=std::min(amin,e);
            const double rr=std::sqrt(lne*(1.0/amin+1.0/g.alpha));   // exact pair threshold (== ForImageOffsets)
            const rvec3_t ci=radials[i]->GetCenter();
            const rvec3_t dig=ci-g.center;
            dcmplx s(0.0);
            for (const auto& n : A.CellsInSphere(rr+norm(dig)))      // exact bound (keep test centred at dig)
            {
                const rvec3_t Roff=A.ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z)));
                const rvec3_t d=dig-Roff;
                if (d.x*d.x+d.y*d.y+d.z*d.z > rr*rr) continue;       // product < eps -> series converged here
                s += phase(n)*GaussOverlapTerm(i, g, g.center+Roff, Lg);
            }
            b[i]=ns[i]*s;
        }
        return b;
    }
    // The FINITE (home-term-only) <chi_i|g>: the molecule/box limit -- no lattice.
    cvec_t MakeOverlap(const Molecule::LatticeSum1E::GaussianFunction& g) const
    {
        const int Lg=GaussDegree(g);
        cvec_t b(size(), dcmplx(0.0));
        for (auto i:indices()) b[i]=ns[i]*GaussOverlapTerm(i, g, g.center, Lg);
        return b;
    }
    //! g's max polynomial degree (Hermite table sizing) and one raw \f$\langle\chi_i|g\,\text{at}\,c\rangle\f$
    //! term -- shared by the periodic and finite \c MakeOverlap(g) forms.
    static int GaussDegree(const Molecule::LatticeSum1E::GaussianFunction& g)
    {
        int Lg=0;
        for (const auto& t : g.terms) Lg=std::max(Lg, t.p.n+t.p.l+t.p.m);
        return Lg;
    }
    double GaussOverlapTerm(size_t i, const Molecule::LatticeSum1E::GaussianFunction& g,
                            const rvec3_t& c, int Lg) const
    {
        GaussianRF gr(g.alpha, c, Lg);
        double sr=0.0;
        for (const auto& t : g.terms)
            sr += t.c * radials[i]->Overlap2C(gr, pols[i], Polarization(t.p));
        return sr;
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

    // GPW collocate/integrate pair-loop threading.  SERIAL BY DEFAULT so the Si bit-anchors stay
    // byte-identical (the threaded path reduces cross-pair sums in a load-dependent order -> a few-ULP
    // drift, the same reason OpenBLAS is pinned to one thread in the harness).  Opt in for the slow NaF
    // production-grid runs with the env knob GPW_OMP_THREADS>1 (read once).  A per-run env knob rather than
    // OMP_NUM_THREADS because the Si anchors and a threaded NaF sweep share one UTMain binary and cannot be
    // separated by a global harness pin -- exactly the NAF_*/GPW_ILLCOND_ECUT env-knob idiom already in use.
    static int PairThreads()
    {
#ifdef QCHEM_OPENMP
        static const int n=[]{ const char* s=std::getenv("GPW_OMP_THREADS"); int v=s?std::atoi(s):1; return v<1?1:v; }();
        return n;
#else
        return 1;
#endif
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

    //! \brief A SHELL: the contiguous run of basis-function indices \f$[begin,end)\f$ sharing one contracted
    //! radial.  \c PGData::Init flattens the basis \c Block by \c Block (one \c GaussianRF + its
    //! polarizations), so a shell IS a contiguous run over which \c radials[i] is the same object -- the
    //! same fact the pointwise sweep's shell hoist uses (PG_Cart::IrrepBasisSet::operator()).
    struct Shell { size_t begin=0, end=0; };
    //! The shell partition of the function index range, built once (geometry-fixed) by scanning \c radials
    //! for runs of the same object.  NOT a pointer-keyed container -- an adjacency scan over a flat vector.
    const std::vector<Shell>& Shells() const
    {
        if (!itsShells.empty() || size()==0) return itsShells;
        const GaussianRF* last=nullptr;
        for (size_t i=0;i<size();i++)
        {
            if (radials[i]!=last) { itsShells.push_back({i,i}); last=radials[i]; }
            itsShells.back().end=i+1;
        }
        return itsShells;
    }

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
        // The enumeration sphere is centred at Roff=0 but the KEEP test at dij, so the exact bound is
        // rr+|dij| (triangle inequality) -- NOT rr+maxCellEdge, whose "|dij| <= cell edge" assumption FAILS
        // for oblique cells / far in-cell pairs (the MnO rhombohedral AFM-II cell: |dij| up to ~2.1x the
        // edge, whole image shells silently missed, indefinite Bloch overlap lambda_min ~ -0.2).
        const double rr=std::sqrt(-std::log(kScreenEps())*(1.0/aMinI+1.0/aMinJ));
        const rvec3_t dij=Ri-Rj;
        for (const auto& n : A.CellsInSphere(rr+norm(dij)))
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
    // The radial's MOST DIFFUSE primitive (smallest exponent = longest reach): the conservative exponent for a
    // magnitude screen (screening on it drops only sub-eps terms).  Used by the 3-centre local-PP screen.
    double MinExponent(size_t i) const
    {
        double a=radials[i]->GetExponents()[0];
        for (double e:radials[i]->GetExponents()) a=std::min(a,e);
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
    double RelCutoffSafety() const {return kRelSafety;}   // exposed via LatticeSum1E (the ladder-completion rung)
    //! \a absRelCutoff = 0: the RELATIVE smooth-field rule (below) -- the density-side calibration, shared by
    //! collocation and the per-iteration integrate-back so the two stay exact adjoints.  \a absRelCutoff > 0:
    //! the ABSOLUTE rule req = absRelCutoff*(alpha_i+alpha_j) (Ha per unit pair exponent -- CP2K's
    //! gaussian_gridlevel REL_CUTOFF; its Ry keyword is 2x this).  The absolute rule bounds the pair's own
    //! spectral tail at its level by e^{-absRelCutoff/2} UNIFORMLY (e^{-15} at 30 Ha), independent of the
    //! FIELD's sharpness -- the property that makes the static local-PP sweep standalone-exact
    //! (doc/GPWPlan.md 0e-PP; the old relCutoffScale multiplier of the relative rule is retired).
    size_t PairLevel(size_t i, size_t j, const std::vector<double>& ecut_L, double absRelCutoff,
                     double fieldSharpness=0.0, double relFieldSharp=-1.0) const
    {
        // ecut_L[0] is the RESOLUTION REFERENCE (the charge-calibrated density grid) -- the relative req is
        // measured against ITS resolution, so appending a finer completion rung (doc/GPWPlan 0b') must not
        // stiffen every pair's requirement.  Selection is order-free: the COARSEST level satisfying req,
        // else the FINEST present (the completion rung when the ladder is complete; the reference grid when
        // it is not).
        // GRID-MATCHING OVERRIDE (doc/GPWPlan §0e; verification instrument): GPW_RELCUTOFF=<Ha> forces the
        // ABSOLUTE rule at the given kappa for EVERY assignment (density side included) -- the CP2K-matching
        // experiment knob.
        static const double kEnvRelCutoff = [](){ const char* s=std::getenv("GPW_RELCUTOFF"); return s ? std::atof(s) : 0.0; }();
        const double kappa = kEnvRelCutoff>0.0 ? kEnvRelCutoff : absRelCutoff;
        // FIELD SHARPNESS in BOTH rules (doc/GPWPlan1.md 4b): the integrand chi_i * V * chi_j is a product of
        // Gaussians, so exponents ADD -- the grid must resolve alpha_i+alpha_j + beta, where beta is the FIELD's
        // own effective exponent, NOT the pair alone.  Without it a DIFFUSE pair (alpha_i+alpha_j tiny) against a
        // SHARP field lands on a coarse level that cannot resolve the field -> a spurious low diffuse eigenvalue.
        //   * ABSOLUTE rule (local PP): beta = fieldSharpness = 1/(2 rloc^2), passed by MakeLocalPP.
        //   * RELATIVE rule (density/Hartree/XC): the KS field's core sharpness follows the valence density,
        //     whose tightest Gaussian product exponent is 2*alpha_max; the LDA XC potential V_xc ~ rho^{1/3}
        //     softens that to ~ (2/3)*alpha_max.  So beta_field = kFieldSharp*MaxExponent() with kFieldSharp=2/3
        //     -- DERIVED, not tuned: it floors every pair's requirement at the resolution the sharp XC core
        //     needs, curing the NaF 3^3-grid diving-ghost the 63f20bd1 coarse ladder exposed.  (GPW_FIELDSHARP
        //     overrides the 2/3 for the self-convergence check.)  0 reproduces the old pair-only relative rule.
        // MAX not ADD: a diffuse pair * sharp field is as sharp as the FIELD alone (the sharp factor dominates
        // the product envelope), so max() LIFTS the diffuse pairs that cannot resolve the field while leaving
        // tight/mid pairs -- already fine enough for their own product -- untouched (no cost blow-up on a
        // diffuse-heavy basis; doc/GPWPlan1.md 4b).  beta_field = kFieldSharp*alpha_max is the KS field's core
        // exponent (V_xc ~ rho^{1/3}); kFieldSharp pinned by rho_lost/N grid-convergence, not wall-clock.
        // relFieldSharp<0 (the default) = the historical beta = kFieldSharp*alpha_max; >=0 = the caller's
        // EXPLICIT beta -- 0 is pair-only routing (RasterFields::HartreeOnly: the raster serves only the
        // smoothing Poisson solve, so a diffuse pair's own bandwidth is the whole requirement).
        static const double kFieldSharp = [](){ const char* s=std::getenv("GPW_FIELDSHARP"); return s?std::atof(s):(2.0/3.0); }();
        const double beta = relFieldSharp>=0.0 ? relFieldSharp : kFieldSharp*MaxExponent();
        const double req = kappa>0.0
            ? kappa*std::max(MaxExponent(i)+MaxExponent(j), fieldSharpness)
            : kRelSafety*ecut_L[0]*std::max(MaxExponent(i)+MaxExponent(j), beta)/(2.0*MaxExponent());
        size_t L=0;
        for (size_t l=1; l<ecut_L.size(); l++) if (ecut_L[l]>ecut_L[L]) L=l;          // fallback: finest present
        bool sat=false;
        for (size_t l=0; l<ecut_L.size(); l++)
            if (ecut_L[l]>=req && (!sat || ecut_L[l]<ecut_L[L])) { L=l; sat=true; }   // coarsest satisfying
        return L;
    }
    //! STATIC-EXTERNAL-FIELD pair->level assignment (LatticeSum1E face; the physics doc there).  The
    //! requirement is the HARMONIC rule 2 lnEps p beta/(p+beta) -- the exact ball that bounds the
    //! truncation SUMMAND e^{-G^2(1/4p+1/4beta)} of a pair p against a static field beta by eps -- floored
    //! by the pair's OWN quadrature requirement (the relative rule, no field floor: the field term already
    //! covers sharp fields, and a static field needs no density-side low-G floor).  Selection = coarsest
    //! satisfying, else finest present (same convention as PairLevel).
    std::vector<size_t> StaticFieldPairLevels(const std::vector<double>& ecut_L,
                                              double beta, double lnEps) const
    {
        assert(beta>0.0 && lnEps>0.0);
        const size_t n=size();
        size_t Lf=0;                                                       // finest present (the fallback)
        for (size_t l=1; l<ecut_L.size(); l++) if (ecut_L[l]>ecut_L[Lf]) Lf=l;
        std::vector<size_t> lv(n*n, Lf);
        for (auto i:indices()) for (auto j:indices(i))
        {
            const double p=MaxExponent(i)+MaxExponent(j);
            const double req=std::max(kRelSafety*ecut_L[0]*p/(2.0*MaxExponent()),   // the pair's own quadrature
                                      2.0*lnEps*p*beta/(p+beta));                   // static-field truncation
            size_t L=Lf;
            bool sat=false;
            for (size_t l=0; l<ecut_L.size(); l++)
                if (ecut_L[l]>=req && (!sat || ecut_L[l]<ecut_L[L])) { L=l; sat=true; }
            lv[i*n+j]=L;
        }
        return lv;
    }
    // Iterate the (i, j@Roff) product's compact exp-tail box on the N-division grid of cell A, visiting each
    // screened-in, MODULO-WRAPPED grid point of chi_i(r) chi_j(r-R_j-Roff).  Shared by
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
    //! \a epsEff: the value tolerance this COLLOCATION box must honour.  Default = kDensityEps (the
    //! collocation accuracy), NOT kScreenEps (the analytic-integral reach) -- the two are decoupled so a
    //! looser GPW_DENSITY_EPS shrinks the collocation boxes for speed WITHOUT loosening S/T/V_local (which
    //! keeps the overlap PSD).  The D-aware path passes max(kDensityEps, kDensityEps/|weight|) -- a
    //! small-weight pair keeps a smaller box (or none: the prefactor early-out below is the whole-term kill).
    //!
    //! \b SHELL-BLOCKED (2026-08-19, doc/GPWPlan1.md "Round 4").  This is the SHELL-pair kernel; the
    //! single-pair \c ForPairBox below is the \f$1\times1\f$ case of it, so there is one box walk in the
    //! codebase.  EVERYTHING here is a SHELL property -- it reads \c radials[i0] / \c radials[j0] only: the
    //! exponents and coefficients, the centres, the product centre \f$P\f$, \a reach, the fractional
    //! half-widths, the ellipsoid pre-screen, the incremental \f$r\f$ walk and the modulo wrap.  A shell's
    //! components differ ONLY in \c pols[i] and \c ns[i].  So the walk runs once per shell pair and each
    //! surviving point hands the caller the two per-component FACTOR arrays
    //! \f$f^I_a=n_{i_0+a}\,\mathrm{pol}_{i_0+a}(d_I)\,\mathrm{rad}_I\f$ and \f$f^J_b\f$ likewise: a
    //! \f$d\times d\f$ shell pair then evaluates 2 contracted radials + 12 polynomials per point instead of
    //! 72 + 72.  The caller forms \f$val=f^I_a f^J_b\f$ for the component pairs it wants and applies its own
    //! \f$|val|<\varepsilon\f$ screen -- the product expression is associated exactly as the unblocked form
    //! was, so values are BIT-IDENTICAL.
    //! \param f called \c f(rasterIndex, fI, fJ) at each point surviving the ellipsoid pre-screen.
    template <class F>
    void ForShellPairBox(size_t i0, size_t nI, size_t j0, size_t nJ, const rvec3_t& Roff,
                         const UnitCell& A, const ivec3_t& N, F&& f, double epsEff=kDensityEps()) const
    {
        const size_t i=i0, j=j0;                             // the shell's representative (radials are shared)
        const rvec_t ei=radials[i]->GetExponents(), gi=radials[i]->GetCoeffs();
        const rvec_t ej=radials[j]->GetExponents(), gj=radials[j]->GetCoeffs();
        const rvec3_t Ri=radials[i]->GetCenter(), Rj=radials[j]->GetCenter()+Roff;   // chi_j imaged to cell Roff
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
        // The per-component factor arrays, allocated ONCE per box (kMaxShell bounds a Cartesian shell:
        // (L+1)(L+2)/2 = 45 at L=8, well past any basis this evaluator sees).
        static constexpr size_t kMaxShell=45;
        assert(nI<=kMaxShell && nJ<=kMaxShell && "ForShellPairBox: shell wider than kMaxShell components");
        double fI[kMaxShell], fJ[kMaxShell];
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
                    // Associated exactly as the unblocked (ni*pols[i](di)*radI)*(nj*pols[j](dj)*radJ) was.
                    for (size_t a=0;a<nI;a++) fI[a]=ns[i0+a]*pols[i0+a](di)*radI;
                    for (size_t b=0;b<nJ;b++) fJ[b]=ns[j0+b]*pols[j0+b](dj)*radJ;
                    const long gx=cx+dx, gy=cy+dy, gz=cz+dz;     // unwrapped grid index (may lie outside [0,N))
                    const long mx=((gx%N.x)+N.x)%N.x, my=((gy%N.y)+N.y)%N.y, mz=((gz%N.z)+N.z)%N.z;  // wrap
                    f((size_t(mx)*N.y+my)*N.z+mz, fI, fJ);
                }
            }
        }
    }
    //! \brief The shell pair's product PREFACTOR exponent at cross-cell offset \a Roff:
    //! \f$\mathrm{pf}=\tfrac{\alpha^I_{\min}\alpha^J_{\min}}{\alpha^I_{\min}+\alpha^J_{\min}}|R_i-R_j-R_{off}|^2\f$,
    //! screen (1) of \c ForShellPairBox.  \f$e^{-\mathrm{pf}}\f$ bounds EVERY value in the box, so a
    //! consumer whose tolerance is \f$\varepsilon\f$ can drop the whole (pair, offset) term when
    //! \f$\mathrm{pf}\ge-\ln\varepsilon\f$ -- without walking anything.
    //!
    //! WHY IT IS EXPOSED (stage B): \c pf is a SHELL property but the D-aware tolerance
    //! \f$\varepsilon_{ij}=\varepsilon/|c_{ij}|\f$ is PER COMPONENT PAIR, so a blocked walk that only
    //! consulted its union tolerance would silently WALK terms the unblocked per-component box killed
    //! outright -- which is what turned the shell hoist into a wash on the D-aware SCF path (measured:
    //! 1.03-1.08x before this, against 2.13x on the uniform-tolerance local-PP sweep).  Pre-filtering the
    //! component pairs on their OWN tolerance restores the kill exactly, and shrinks the union box to the
    //! survivors as a bonus.
    double PairPrefactorExp(size_t i0, size_t j0, const rvec3_t& Roff) const
    {
        const rvec_t ei=radials[i0]->GetExponents(), ej=radials[j0]->GetExponents();
        double aMinI=ei[0]; for (double e:ei) aMinI=std::min(aMinI,e);
        double aMinJ=ej[0]; for (double e:ej) aMinJ=std::min(aMinJ,e);
        const rvec3_t dij=radials[i0]->GetCenter()-(radials[j0]->GetCenter()+Roff);
        return (aMinI*aMinJ/(aMinI+aMinJ))*(dij.x*dij.x+dij.y*dij.y+dij.z*dij.z);
    }
    //! The SINGLE-PAIR box (the \f$1\times1\f$ case of \c ForShellPairBox): calls
    //! \c f(rasterIndex, \f$\chi_i(r)\chi_j(r-R_j-R_{off})\f$) at each screened-in, wrapped grid point.
    //! Kept as the name every single-pair consumer speaks; the walk itself is not duplicated.
    template <class F>
    void ForPairBox(size_t i, size_t j, const rvec3_t& Roff, const UnitCell& A, const ivec3_t& N, F&& f,
                    double epsEff=kDensityEps()) const
    {
        ForShellPairBox(i,1,j,1,Roff,A,N,
                        [&](size_t idx, const double* fI, const double* fJ)
                        {
                            const double val=fI[0]*fJ[0];
                            if (std::fabs(val)<epsEff) return;
                            f(idx,val);
                        }, epsEff);
    }
    // --- STREAM CACHE for the analytic collocation / integrate-back ----------------------------------------
    // The per-(pair, cross-cell offset) box streams -- wrapped grid RUNS + analytic pair values -- are pure
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
    //!
    //! \b RUN-LENGTH \b GEOMETRY (2026-08-19).  The raster indices are not stored per point: \c ForPairBox
    //! walks the box's innermost axis (\f$z\f$) in unit steps, so the wrapped indices it emits come in long
    //! CONTIGUOUS runs (one per \f$(dx,dy)\f$ column, split only where the modulo wrap or the value screen
    //! interrupts it).  A stream therefore stores \c runBase/\c runLen -- (first index, length) per run --
    //! and the values alone; the replay walks \c val sequentially and the destination grid contiguously
    //! (\c dst[t] \c += \c cw*v[k]).
    //! WHY: the replay is DRAM-BANDWIDTH bound, not arithmetic bound (perf 2026-08-19 on the MnO magnetic
    //! cell: 49% of \c CollocateDensity's self time sat on the two sequential loads -- the index stream and
    //! the value stream -- against ~1% on the arithmetic and 0% on the Bloch-phase contraction).  Per-point
    //! indices doubled the fp32 tier's traffic (4 B value + 4 B index) and added a third of the fp64 tier's;
    //! a run head costs 8 B once per RUN instead, so the encoding pays whenever the mean run exceeds 2
    //! points.  Measured \c meanRun: 10.2 on the MnO magnetic cell (647 M points in 63.7 M runs -- tiers
    //! 12 -> 8.8 and 8 -> 4.8 B/pt, stream RAM 5.78 -> 3.70 GB) and 4.2 on the small Si crystal gate.
    //! Replay order, values and screens are UNCHANGED, so every replay stays bit-identical.
    struct PairOffsetStream
    {
        ivec3_t n;
        std::vector<unsigned> runBase;   //!< first raster index of each contiguous run
        std::vector<unsigned> runLen;    //!< its length; \f$\sum\f$ runLen == the stream's point count
        std::vector<double> val;         //!< tier 1 values, run-major (exactly one of val/val32 is filled)
        std::vector<float>  val32;       //!< tier 2 values
        double maxv=0.0;
        //! The stream's point count (the budget currency) -- the filled value tier's length.
        size_t Points() const {return val.empty() ? val32.size() : val.size();}
        //! Append one screened-in point: extend the open run when \a idx continues it, else open a new one.
        void Append(size_t idx, double v, bool f32)
        {
            if (!runLen.empty() && size_t(runBase.back())+runLen.back()==idx) runLen.back()++;
            else { runBase.push_back(unsigned(idx)); runLen.push_back(1); }
            if (f32) val32.push_back(float(v)); else val.push_back(v);
            maxv=std::max(maxv,std::fabs(v));
        }
    };
    struct PairStreams     { size_t level=0; bool cached=false; std::vector<PairOffsetStream> offsets; };
    struct StreamCache
    {
        std::vector<ivec3_t> N_L; std::vector<double> ecut_L;   // the ladder shape (snapshot key)
        double absRelCutoff=0.0;                                //   + the assignment rule key (0=relative, >0=absolute Ha/exponent)
        double relFieldSharp=-1.0;                              //   + the relative rule's beta floor (part of the level-assignment identity)
        std::vector<PairStreams> pairs;                         // indexed [i*n+j], j>=i only
        size_t droppedPts=0;                                    // points that fell past BOTH tiers at build time
        size_t availAtBuild64=0, availAtBuild32=0;              // the budget headroom this build was offered --
                                                                //   an INCOMPLETE cache rebuilds when headroom grows
                                                                //   (a resident shape was released; doc/GPWPlan 0.5(b))
        bool SameShape(const std::vector<ivec3_t>& N, const std::vector<double>& e) const
        {
            if (ecut_L!=e || N_L.size()!=N.size()) return false;
            for (size_t l=0;l<N.size();l++)
                if (N_L[l].x!=N[l].x || N_L[l].y!=N[l].y || N_L[l].z!=N[l].z) return false;
            return true;
        }
    };
    static constexpr size_t kMaxStreamCaches=4;                 // ladder + K=1 + unit-gate shapes; guards churn
    //! Global stream-cache budgets in POINTS, TWO TIERS.  A diffuse basis in a small cell keeps hundreds of
    //! screened offsets per pair -- caching them ALL blew past physical RAM on NaF (uncapped: 27 GB virt on a
    //! 16 GB box; measured demand 952M pts).  Tier 1 (fp64, ~8.8 B/pt since the run-length encoding,
    //! ~1.3 GB): replay is BIT-IDENTICAL to
    //! on-the-fly evaluation; sized so the Si SR ladder caches completely (demand 104.9M) -- every Si anchor
    //! and machine-precision kernel gate lives entirely in this tier.  Tier 2 (fp32 values, ~4.8 B/pt, ~4.1 GB):
    //! the overflow pairs store float values instead of falling to on-the-fly -- replay noise is ~6e-8
    //! RELATIVE per value (NaF anchors are 1e-3-scale; Tr(DS) charge is analytic, unaffected), and the
    //! collocate/integrate ADJOINT stays machine-exact because both directions replay the SAME stream.
    //! Pairs beyond BOTH budgets fall back to on-the-fly evaluation (correct, slower).
    static constexpr size_t kStreamBudgetPts   =150'000'000;   // tier 1: fp64, bit-identical replay
    //! Tier 2 sized so NaF's MEASURED demand (952M total: 150M fp64 + 802M fp32) caches COMPLETELY -- at
    //! 700M its last 314 (small) pairs re-evaluated on the fly every sweep.  Peak stream RAM ~ 8.6 GB.
    static constexpr size_t kStreamBudgetPtsF32=850'000'000;   // tier 2: fp32 values (the NaF coverage lever)
    //! RUNTIME budget overrides (POINTS; 0/unset = the compile-time defaults above).  A MEMORY-SAFETY valve
    //! for boxes where the defaults' ~8.6 GB peak cannot fit beside the rest of a run (measured 2026-07-22:
    //! the FULL VALENCE_LOWQ_SR basis -- alpha_min=0.05, reach ~21.5 au, hundreds of offsets/pair -- drove
    //! UTMain to 12.6 GB on a 14 GB box and the kernel OOM-killed the desktop).  Dropped pairs fall back to
    //! on-the-fly evaluation: correct, slower -- never a physics knob.
    static size_t BudgetPts()    { const char* e=std::getenv("GPW_STREAM_BUDGET_PTS");     return e?size_t(std::atoll(e)):kStreamBudgetPts; }
    static size_t BudgetPtsF32() { const char* e=std::getenv("GPW_STREAM_BUDGET_PTS_F32"); return e?size_t(std::atoll(e)):kStreamBudgetPtsF32; }
    //! The budget headroom a NEW build would be offered: the global tier budgets minus the points every
    //! cache EXCEPT \a skip already holds (budgets are GLOBAL across cache shapes).
    std::pair<size_t,size_t> StreamBudgetHeadroom(const StreamCache* skip) const
    {
        size_t a64=BudgetPts(), a32=BudgetPtsF32();
        for (const auto& cc : itsStreamCaches)
        {
            if (&cc==skip) continue;
            for (const auto& ps : cc.pairs)
                for (const auto& st : ps.offsets)
                    if (!st.val.empty()) a64 -= std::min(a64, st.Points());
                    else                 a32 -= std::min(a32, st.Points());
        }
        return {a64,a32};
    }
    //! doc/GPWPlan.md 0.5(b): drop the caches for one ladder shape (any assignment rule), refunding the
    //! global budget.  Realises Molecule::LatticeSum1E::ReleaseStreams (forwarded by the host IBS); called
    //! by a dying grid client (GPW_Evaluator dtor) so a finished stage's streams stop squatting on the
    //! budget.  An INCOMPLETE surviving cache picks the refund up via the self-heal rebuild in EnsureStreams.
    void ReleaseStreams(const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L) const
    {
        for (auto ci=itsStreamCaches.begin(); ci!=itsStreamCaches.end(); )
            if (ci->SameShape(N_L,ecut_L)) ci=itsStreamCaches.erase(ci);
            else ++ci;
    }
    // --- T3 route (b) STREAM FOLD (doc/SymmetryUpgradePlan.md §6b) -----------------------------------------
    // Fold the (pair i<=j, offset n) collocation terms under the imposed crystal ops {W|tau}.  Pure
    // geometry-fixed bookkeeping, cached like the stream caches: per accepted op, the basis map
    // i -> (i', s_i, L_i) (partner function + Cartesian-monomial sign + integer cell offset of the image
    // centre); per canonical pair slot, either REPRESENTATIVE (with the pair-stabilizer's action on its
    // offset list) or IMAGE (edge to its rep: sigma = s_i s_j, Hermitian flip).  A triple maps as
    //     g.(i,j,n) = (i', j', W n + L_j - L_i),   sign sigma = s_i s_j,
    // canonicalized through the Hermitian twin (i,j,n)~(j,i,-n); on the DIAGONAL the twin acts within the
    // pair ((i,i,n) and (i,i,-n) have the SAME wrapped product field -- an integer-cell translate), carried
    // as an extra offset involution.  For a group-symmetric D (D_{i'j'} = sigma D_ij -- what §3 imposition
    // asserts) the signs cancel between weight and stream, so REDUCED replay scatters each representative
    // triple with its plain orbit multiplicity (pairMult x offset within-multiplicity) and the caller's
    // dense group-average (SymmetrizeGMap / SymmetrizeRaster -- exact at ANY raster N: integer W = exact
    // voxel map, tau rides the FFT shift) reproduces the full collocation.  NO sparse stream is ever
    // index-permuted, hence NO §5 raster-commensurability precondition (the §6b verdict).  Integrate-back
    // gathers representative pairs only (within-multiplicity weights) and fills partners by the
    // representation transform h_{i'j'} = sigma h_ij -- exact for a group-symmetric V, and the exact
    // adjoint of the reduced scatter (symmetrize V first; §6b item 5).  A pair fixed by an op with
    // sigma = -1 carries D_ij = 0 under the imposed group: DEAD (skipped; h = 0 is its projected value).
    // GENERAL k (T3.4): the fold takes the block's kFrac and folds under the LITTLE GROUP of k only
    // (guarded here; SpaceGroup::LittleGroupDirectOps is the caller's source).  For little-group ops
    // e^{2πik·Wn} == e^{2πik·n} (mod 1), so the replay multiplicities stay plain integers; k enters only
    // the edge factor zeta = sigma e^{2πik·(L_j-L_i)} -- the dead rule (non-flip self-edge with zeta != 1)
    // and the h image fill h_{i'j'} = zeta h_ij.  Flip (Hermitian-twin) self-edges at k != Γ pin only the
    // PHASE of D_ij (conj(D) = zeta D): never dead, conservatively not folded under.
    struct OpMap     { Matrix3D<double> W; std::vector<unsigned> to; std::vector<signed char> sgn; std::vector<ivec3_t> L; };
    struct StabEntry { Matrix3D<double> W; ivec3_t d; bool neg; };            // offset action n -> ±(W n + d)
    //! \a zeta is the edge factor \f$\zeta=\sigma\,e^{2\pi i\,k\cdot(L_j-L_i)}\f$ (§6b/T3.4): the imposed
    //! constraint reads \f$D^k_{i'j'}=\zeta D^k_{ij}\f$ and the image fill \f$h_{i'j'}=\zeta h_{ij}\f$.
    //! At \f$\Gamma\f$ it is the real monomial sign \f$\sigma=\pm1\f$.
    struct PairEdge  { int rep=-1; dcmplx zeta=dcmplx(1.0); bool flip=false; bool dead=false; unsigned pairMult=1; };
    struct StreamFold
    {
        std::vector<Symmetry::Lattice_3D::DirectOp> srcOps; //!< the caller's op set (idempotence key -- k-blocks re-inject)
        rvec3_t                                   kF;       //!< the fractional k the fold was built for (idempotence key)
        std::vector<OpMap>                        maps;     //!< accepted ops' basis actions
        std::vector<PairEdge>                     pairs;    //!< [i*n+j] (j>=i slots); rep==own slot marks a representative
        std::vector<std::vector<StabEntry>>       stab;     //!< rep slots: the pair-stabilizer's offset action (identity included)
        std::vector<std::map<long long,unsigned>> offMult;  //!< rep slots: offset key -> within-pair orbit multiplicity (0 = member)
    };
    static long long OffKey(const ivec3_t& v)
    {   return ((long long)(v.x+(1<<20))<<42) | ((long long)(v.y+(1<<20))<<21) | (long long)(v.z+(1<<20)); }
    static ivec3_t IntApply(const Matrix3D<double>& W, const ivec3_t& v)     // exact integer action (W integer-valued)
    {
        return ivec3_t(int(std::lround(W(1,1)*v.x+W(1,2)*v.y+W(1,3)*v.z)),
                       int(std::lround(W(2,1)*v.x+W(2,2)*v.y+W(2,3)*v.z)),
                       int(std::lround(W(3,1)*v.x+W(3,2)*v.y+W(3,3)*v.z)));
    }
    //! Decompose the op's CARTESIAN action \f$R=A\,W\,A^{-1}\f$ as a signed axis permutation
    //! \f$(R^\top u)_a = s_a\,u_{q_a}\f$.  False when \f$R\f$ genuinely mixes axes (non-cubic lattices) --
    //! the op is then dropped from the fold (folding merely reduced; general shell-mixing is a later
    //! increment, doc/SymmetryUpgradePlan.md §6b item 2).
    bool CartSignedPerm(const UnitCell& A, const Matrix3D<double>& W, int q[3], int s[3]) const
    {
        for (int a=0;a<3;a++)
        {
            const rvec3_t e(a==0?1.0:0.0, a==1?1.0:0.0, a==2?1.0:0.0);
            const rvec3_t col=A.ToCartesian(W*A.ToFractional(e));            // column a of R
            const double c[3]={col.x,col.y,col.z};
            int nz=-1;
            for (int b=0;b<3;b++)
                if (std::fabs(c[b])>1e-9)
                {
                    if (nz>=0 || std::fabs(std::fabs(c[b])-1.0)>1e-9) return false;
                    nz=b;
                }
            if (nz<0) return false;
            q[a]=nz; s[a]=c[nz]>0?1:-1;
        }
        return q[0]!=q[1] && q[1]!=q[2] && q[0]!=q[2];
    }
    //! Realises \c Molecule::LatticeSum1E::SetStreamSymmetryOps (forwarded by the host IBS): build (or
    //! clear, \a ops empty) the stream fold for a block at fractional crystal momentum \a kFrac
    //! (default \f$\Gamma\f$).  Returns the number of ops actually used.  Const like the stream caches --
    //! the fold is derived, geometry-fixed bookkeeping.
    //! T3.4: ops must lie in the LITTLE GROUP of \a kFrac (\c SpaceGroup::LittleGroupDirectOps); a
    //! non-little-group op is DROPPED here (belt + braces -- it relates different k-blocks, not terms
    //! within this one).  For little-group ops the replay weights keep their plain integer
    //! multiplicities (the \f$e^{2\pi ik\cdot Wn}\equiv e^{2\pi ik\cdot n}\f$ cancellation); \a kFrac
    //! enters only the edge factors \f$\zeta\f$ (dead rule + h image fill).
    size_t SetStreamSymmetryOps(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops, const UnitCell& A,
                                const rvec3_t& kFrac=rvec3_t(0,0,0)) const
    {
        // IDEMPOTENT: sibling k-blocks re-inject the same set (the evaluator is shared); an identical set
        // keeps the existing fold AND the stream caches built under it.  Any CHANGE (including clearing)
        // invalidates the caches -- a REDUCED-built cache must never serve a free-run replay (§6b/T3.2).
        auto sameOps=[&]()->bool
        {
            if (!itsStreamFold || itsStreamFold->srcOps.size()!=ops.size()) return false;
            const rvec3_t& k0=itsStreamFold->kF;
            if (k0.x!=kFrac.x || k0.y!=kFrac.y || k0.z!=kFrac.z) return false;
            for (size_t o=0;o<ops.size();o++)
            {
                const auto &a=itsStreamFold->srcOps[o], &b=ops[o];
                for (size_t r=1;r<=3;r++) for (size_t c=1;c<=3;c++) if (a.W(r,c)!=b.W(r,c)) return false;
                if (a.tau.x!=b.tau.x || a.tau.y!=b.tau.y || a.tau.z!=b.tau.z) return false;
            }
            return true;
        };
        if (sameOps()) return itsStreamFold->maps.size();
        if (itsStreamFold || !ops.empty()) itsStreamCaches.clear();
        itsStreamFold.reset();
        if (ops.empty()) return 0;
        auto sf=std::make_unique<StreamFold>();
        sf->srcOps=ops;
        sf->kF=kFrac;
        const bool kZero = kFrac.x==0.0 && kFrac.y==0.0 && kFrac.z==0.0;
        const double twoPi=8.0*std::atan(1.0);
        const size_t nn=size();
        std::vector<rvec3_t> fc(nn);                                         // fractional centres
        for (auto i:indices()) fc[i]=A.ToFractional(radials[i]->GetCenter());
        auto sameRadial=[&](size_t a, size_t b)->bool
        {
            const rvec_t &ea=radials[a]->GetExponents(), &eb=radials[b]->GetExponents();
            const rvec_t &ga=radials[a]->GetCoeffs(),    &gb=radials[b]->GetCoeffs();
            if (ea.size()!=eb.size() || ga.size()!=gb.size()) return false;
            for (size_t p=0;p<ea.size();p++) if (std::fabs(ea[p]-eb[p])>1e-10*(1.0+std::fabs(ea[p]))) return false;
            for (size_t p=0;p<ga.size();p++) if (std::fabs(ga[p]-gb[p])>1e-10*(1.0+std::fabs(ga[p]))) return false;
            return true;
        };
        for (const auto& op : ops)
        {
            if (!kZero)                                                      // little-group guard (T3.4):
            {                                                                //   (W^{-1})^T k must == k (mod 1)
                const rvec3_t uk=Transpose(Invert(op.W))*kFrac;
                const rvec3_t dk(uk.x-kFrac.x, uk.y-kFrac.y, uk.z-kFrac.z);
                if (std::fabs(dk.x-std::round(dk.x))>1e-9 || std::fabs(dk.y-std::round(dk.y))>1e-9
                                                          || std::fabs(dk.z-std::round(dk.z))>1e-9) continue;
            }
            int q[3], sax[3];
            if (!CartSignedPerm(A,op.W,q,sax)) continue;                     // shell-mixing op: dropped
            OpMap m; m.W=op.W; m.to.resize(nn); m.sgn.resize(nn); m.L.resize(nn);
            bool ok=true;
            for (size_t i=0;i<nn && ok;i++)
            {
                const rvec3_t ft=op.W*fc[i]+op.tau;                          // image centre (fractional)
                Polarization pm;                                             // permuted monomial + sign
                int sg=1;
                for (int a=0;a<3;a++)
                {
                    pm[q[a]]=pols[i][a];
                    if (sax[a]<0 && (pols[i][a]&1)) sg=-sg;
                }
                ok=false;
                for (size_t j=0;j<nn;j++)
                {
                    if (pols[j]!=pm) continue;
                    const rvec3_t d=ft-fc[j];
                    const double Lx=std::round(d.x), Ly=std::round(d.y), Lz=std::round(d.z);
                    if (std::fabs(d.x-Lx)>1e-6 || std::fabs(d.y-Ly)>1e-6 || std::fabs(d.z-Lz)>1e-6) continue;
                    if (!sameRadial(i,j)) continue;
                    m.to[i]=unsigned(j); m.sgn[i]=(signed char)sg; m.L[i]=ivec3_t(int(Lx),int(Ly),int(Lz));
                    ok=true; break;
                }
            }
            if (ok) sf->maps.push_back(std::move(m));                       // op maps the basis onto itself
        }
        if (sf->maps.empty()) return 0;
        // Pair fold: process canonical slots ascending -- the first unassigned slot of an orbit IS its rep
        // (orbits are closed under the accepted op set, one application each covers the orbit).
        sf->pairs.assign(nn*nn,{});
        sf->stab.assign(nn*nn,{});
        for (auto i:indices()) for (auto j:indices(i))
        {
            const size_t s0=i*nn+j;
            if (sf->pairs[s0].rep>=0) continue;
            sf->pairs[s0]={int(s0),dcmplx(1.0),false,false,1};
            bool dead=false;
            std::vector<size_t> members{s0};
            for (const auto& m : sf->maps)
            {
                unsigned x=m.to[i], y=m.to[j];
                const int sg=int(m.sgn[i])*int(m.sgn[j]);
                const ivec3_t d(m.L[j].x-m.L[i].x, m.L[j].y-m.L[i].y, m.L[j].z-m.L[i].z);
                // The edge factor zeta = sigma e^{2πi k.(L_j-L_i)} (§6b/T3.4): D^k_{i'j'} = zeta D^k_ij.
                const dcmplx zeta = kZero ? dcmplx(double(sg))
                                          : double(sg)*std::polar(1.0, twoPi*(kFrac.x*d.x+kFrac.y*d.y+kFrac.z*d.z));
                bool fl=false;
                if (x>y) { std::swap(x,y); fl=true; }
                const size_t s1=size_t(x)*nn+y;
                if (s1==s0)
                {
                    // SELF-EDGES.  Non-flip: D_ij = zeta D_ij, so zeta != 1 kills the pair (the Γ sigma=-1
                    // rule generalized); zeta == 1 folds its offsets.  Flip: the constraint is
                    // conj(D_ij) = zeta D_ij -- a PHASE pin, not deadness, except at Γ where D is real
                    // (then zeta=-1 does kill).  At k != Γ flip self-edges are conservatively SKIPPED
                    // (no offset folding under them -- folding merely reduced, never wrong).
                    if (!fl)
                    {
                        if (std::abs(zeta-dcmplx(1.0))>1e-9) dead=true;
                        else sf->stab[s0].push_back({m.W,d,fl});
                    }
                    else if (kZero)
                    {
                        if      (sg<0) dead=true;                           // Γ, real D: conj==id, zeta=-1 => D=0
                        else           sf->stab[s0].push_back({m.W,d,fl});
                    }
                }
                else if (sf->pairs[s1].rep<0)
                {
                    sf->pairs[s1]={int(s0),zeta,fl,false,1};
                    members.push_back(s1);
                }
                else if (sf->pairs[s1].flip==fl)
                {
                    // Two same-flavour edges to one image must agree; disagreement forces D = 0 (dead).
                    if (std::abs(sf->pairs[s1].zeta-zeta)>1e-9) dead=true;
                }
                else if (kZero && std::abs(sf->pairs[s1].zeta-zeta)>1e-9) dead=true;
                // (k != Γ, flip vs non-flip edges to the same image: a phase pin on D, never deadness.)
            }
            if (i==j)                                                       // diagonal Hermitian twin: n ~ -n
            {
                const size_t k0=sf->stab[s0].size();
                for (size_t e=0;e<k0;e++)
                {   StabEntry t=sf->stab[s0][e]; t.neg=!t.neg; sf->stab[s0].push_back(t); }
            }
            sf->pairs[s0].pairMult=unsigned(members.size());
            if (dead) for (size_t s : members) sf->pairs[s].dead=true;
        }
        // Offset fold per representative pair: the first-ENUMERATED member of each within-pair orbit is its
        // rep (deterministic -- build, replay and the on-the-fly path share the ForImageOffsets order).
        sf->offMult.assign(nn*nn,{});
        for (auto i:indices()) for (auto j:indices(i))
        {
            const size_t s0=i*nn+j;
            const PairEdge& pe=sf->pairs[s0];
            if (pe.dead || pe.rep!=int(s0)) continue;
            auto& mm=sf->offMult[s0];
            const auto& st=sf->stab[s0];
            ForImageOffsets(i,j,A,[&](const ivec3_t& nOff, const rvec3_t&)
            {
                if (mm.count(OffKey(nOff))) return;                         // folded into an earlier orbit
                std::vector<ivec3_t> orb{nOff};
                for (const auto& e : st)
                {
                    ivec3_t img=IntApply(e.W,nOff);
                    img=ivec3_t(img.x+e.d.x, img.y+e.d.y, img.z+e.d.z);
                    if (e.neg) img=ivec3_t(-img.x,-img.y,-img.z);
                    bool seen=false;
                    for (const auto& o : orb) if (o.x==img.x && o.y==img.y && o.z==img.z) { seen=true; break; }
                    if (!seen) orb.push_back(img);
                }
                mm[OffKey(nOff)]=unsigned(orb.size());
                for (size_t k=1;k<orb.size();k++) mm.emplace(OffKey(orb[k]),0u);
            });
        }
        const size_t used=sf->maps.size();
        itsStreamFold=std::move(sf);
        return used;
    }
    //! Realises \c Molecule::LatticeSum1E::StreamFoldOrder (the fold-state cache-key input).
    size_t StreamFoldOrder() const {return itsStreamFold ? itsStreamFold->maps.size() : 0;}

    //! PROJECT, DON'T SAMPLE (T3.5, doc/SymmetryUpgradePlan.md §6b): the ORBIT AVERAGE of \a D over each
    //! pair orbit, indexed by canonical slot \f$[i\,n+j]\f$ and filled at REPRESENTATIVE slots only (the
    //! only ones the reduced replay reads).  Empty when no fold is armed.
    //!
    //! The reduced replay reads ONE D element per orbit.  Reading the representative's own value SAMPLES
    //! the orbit, which is a projection only if \f$D\f$ already obeys \f$D_{i'j'}=\zeta D_{ij}\f$ -- i.e.
    //! it asserts a symmetric iterate, which a DEGENERATE OPEN SHELL at integer fill breaks permanently
    //! (the 2026-08-03 default-on retraction).  Averaging instead makes the reduced path read
    //! \f$\bar D=P D\f$, and because collocation is LINEAR and EQUIVARIANT
    //! (\f$\rho[g\cdot D]=g\cdot\rho[D]\f$) the caller's dense star-average then gives
    //! \f$P\,\rho_{\rm red}[D]=\rho[P D]=P\,\rho_{\rm full}[D]\f$ -- EXACTLY the unfolded imposed run's
    //! density, for ANY iterate.  So the fold stops imposing more than \c imposeSymmetry already does and
    //! becomes a pure cost decision.  A DEAD pair projects to zero (skipped: its \f$\bar D\f$ IS zero).
    //! Edge algebra: a member slot \f$s\f$ carries \f$D_s=\zeta D_{\rm rep}\f$, or
    //! \f$D_s=\overline{\zeta D_{\rm rep}}\f$ across the Hermitian twin (\c flip), so it contributes
    //! \f$\bar\zeta D_s\f$ resp. \f$\bar\zeta\,\overline{D_s}\f$ to the representative's mean.
    std::vector<dcmplx> FoldProjectedD(const chmat_t& D) const
    {
        const StreamFold* sf=itsStreamFold.get();
        if (!sf) return {};
        const size_t nn=size();
        std::vector<dcmplx> Dp(nn*nn, dcmplx(0.0));
        for (auto i:indices()) for (auto j:indices(i))
        {
            const PairEdge& pe=sf->pairs[i*nn+j];
            if (pe.dead) continue;
            const dcmplx d(D(i,j));
            Dp[size_t(pe.rep)] += std::conj(pe.zeta)*(pe.flip ? std::conj(d) : d);
        }
        for (auto i:indices()) for (auto j:indices(i))
        {
            const size_t s=i*nn+j;
            const PairEdge& pe=sf->pairs[s];
            if (pe.dead || pe.rep!=int(s) || pe.pairMult<2) continue;
            Dp[s]/=double(pe.pairMult);
        }
        return Dp;
    }

    //! The SCREEN's orbit-consistent value: \f$\max_s|{\rm scr}_s|\f$ over each pair orbit, at representative
    //! slots (empty when no fold is armed).  NOT \c FoldProjectedD -- the integrate-back's \a screenD is the
    //! collocation memo's \c Dscr, a matrix of MAGNITUDES (the elementwise max of \f$|D|\f$ over the spin
    //! channels), not a density matrix.  Signed averaging is meaningless on it and actively harmful: an
    //! orbit whose edges carry mixed \f$\sigma=\pm1\f$ cancels to ~0, the D-aware screen then kills terms
    //! that carry real weight, and H comes out wrong (measured 2026-08-19 on the imposed O2-in-box triplet:
    //! a 2.3 Ha variational collapse).  A screen must be orbit-INVARIANT so every member is truncated the
    //! same way, and the only safe invariant reduction of a magnitude is the MAXIMUM -- it can only keep
    //! terms, never drop one a member would have kept (the no-cut discipline).
    std::vector<double> FoldScreenMax(const chmat_t& scr) const
    {
        const StreamFold* sf=itsStreamFold.get();
        if (!sf) return {};
        const size_t nn=size();
        std::vector<double> m(nn*nn, 0.0);
        for (auto i:indices()) for (auto j:indices(i))
        {
            const PairEdge& pe=sf->pairs[i*nn+j];
            if (pe.dead) continue;
            double& dst=m[size_t(pe.rep)];
            dst=std::max(dst, std::abs(dcmplx(scr(i,j))));
        }
        return m;
    }

    const StreamCache& EnsureStreams(const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                     const std::vector<double>& ecut_L, double absRelCutoff,
                                     double relFieldSharp=-1.0) const
    {
        for (auto ci=itsStreamCaches.begin(); ci!=itsStreamCaches.end(); ++ci)
        {
            const StreamCache& c=*ci;
            if (c.absRelCutoff!=absRelCutoff || c.relFieldSharp!=relFieldSharp || !c.SameShape(N_L,ecut_L)) continue;
            if (c.droppedPts==0) return c;                        // complete: replay is bit-stable, never rebuilt
            // SELF-HEAL (doc/GPWPlan 0.5(b)): this cache dropped pairs because OTHER resident caches held the
            // budget when it was built (the grid-continuation seed builds the FINE shape while the coarse
            // stage still squats).  If a release has since grown the headroom, rebuild -- the dropped pairs
            // re-evaluate EVERY iteration, so one rebuild sweep is always the cheaper side of the trade.
            auto [a64,a32]=StreamBudgetHeadroom(&c);
            if (a64<=c.availAtBuild64 && a32<=c.availAtBuild32) return c;
            itsStreamCaches.erase(ci);
            break;                                                // fall through to a fresh build
        }
        if (itsStreamCaches.size()>=kMaxStreamCaches) itsStreamCaches.clear();   // unexpected shape churn
        qchem::report::Timed timer("setup: collocation stream build");           // THE diffuse-basis setup cost
        itsStreamCaches.emplace_back();
        StreamCache& c=itsStreamCaches.back();
        c.N_L=N_L; c.ecut_L=ecut_L; c.absRelCutoff=absRelCutoff; c.relFieldSharp=relFieldSharp;
        const size_t n=size();
        c.pairs.resize(n*n);
        auto [budget64,budget32]=StreamBudgetHeadroom(&c);        // (empty c contributes nothing)
        c.availAtBuild64=budget64; c.availAtBuild32=budget32;
        size_t nPairs=0, nCached64=0, nCached32=0, pts64=0, pts32=0, ptsDropped=0, nOffs=0;
        // T3 route (b) REDUCED BUILD (§6b/T3.2): under an armed fold, only orbit-REPRESENTATIVE pairs and
        // representative offsets are built (the replay never touches the rest) -- demand and build time
        // fall by ~|orbit|, which also lifts pairs INTO the fp64 tier (the §6b measured accuracy note).
        // SetStreamSymmetryOps clears the caches on any fold change, so a reduced cache never serves an
        // unfolded replay (and vice versa).
        const StreamFold* bsf=itsStreamFold.get();
        auto pairSkip=[&](size_t i, size_t j)->bool
        {   return bsf && (bsf->pairs[i*n+j].dead || bsf->pairs[i*n+j].rep!=int(i*n+j)); };
        auto offSkip=[&](size_t i, size_t j, const ivec3_t& nn)->bool
        {
            if (!bsf) return false;
            const auto& mm=bsf->offMult[i*n+j];
            const auto it=mm.find(OffKey(nn));
            return it==mm.end() || it->second==0;
        };
        // ---- SHELL-BLOCKED TWO-PASS BUILD (2026-08-19, doc/GPWPlan1.md "Round 4") ------------------------
        // Everything the box walk needs -- the pair->level assignment, the screened offset list, the box
        // centre/reach/half-widths, the ellipsoid pre-screen, the incremental r walk and the modulo wrap --
        // is a SHELL property (it reads radials[i] only), and a shell's components share one contracted
        // radial.  So the walk runs ONCE PER SHELL PAIR and every component pair takes its value from the
        // two per-point FACTOR arrays ForShellPairBox hands back: a d x d shell pair evaluates 2 contracted
        // radials + 12 polynomials per point instead of 72 + 72.  Values are bit-identical (same expression,
        // same association), so replay is unchanged.
        //
        // WHY TWO PASSES, ALWAYS (this replaced a serial one-pass build beside a threaded two-pass one).
        // The TIERING ORDER is part of the result: each pair's fp64/fp32/drop tier consumes the GLOBAL
        // budgets in ROW-MAJOR COMPONENT-pair order, so a shell-major single pass would tier differently.
        // Splitting into (1) count, (2) a serial row-major budget walk on the counts, (3) build only the
        // tiered pairs keeps the budget walk in its original order while both EVALUATION passes go
        // shell-major -- and it retires the one-pass path's transient-bound abort and its streaming fp32
        // demotion outright, since a two-pass build knows every pair's size before it materialises anything.
        // Cost: ~2x box evaluations (the counting pass IS an evaluation pass -- the |val|>=eps screen needs
        // the values), repaid many times over by the shell hoist.  Threading (opt-in, GPW_OMP_THREADS)
        // partitions both passes over SHELL pairs; a shell pair owns a DISJOINT set of component pairs, so
        // there is no reduction anywhere and serial and threaded builds stay bit-identical.
        const std::vector<Shell>& shells=Shells();
        std::vector<std::pair<size_t,size_t>> sprs;               // shell pairs (a<=b)
        for (size_t a=0;a<shells.size();a++) for (size_t b=a;b<shells.size();b++) sprs.push_back({a,b});
        std::vector<std::pair<size_t,size_t>> prs;                // component pairs, ROW-MAJOR upper triangle
        for (auto i:indices()) for (auto j:indices(i)) prs.push_back({i,j});
        const double eps=kDensityEps();                           // the build's UNIFORM value screen
        std::vector<size_t> ptsAll(n*n,0);                        // per component pair (disjoint across shells)
        std::vector<char>   tierBuild(n*n,0), tierF32(n*n,0);     // filled by the budget walk between the passes
        std::exception_ptr  firstEx;                              // throw containment -- see CollocateDensity

        // The component pairs a shell pair owns (upper triangle), minus the fold's image/dead ones and --
        // when \a tiered -- minus everything the budget walk did not tier.
        auto ownedPairs=[&](size_t sp, bool tiered, std::vector<std::pair<size_t,size_t>>& out)
        {
            const Shell& si=shells[sprs[sp].first]; const Shell& sj=shells[sprs[sp].second];
            out.clear();
            for (size_t i=si.begin;i<si.end;i++)
                for (size_t j=std::max(i,sj.begin);j<sj.end;j++)      // b>a => every j exceeds every i
                    if (!pairSkip(i,j) && (!tiered || tierBuild[i*n+j])) out.push_back({i,j});
        };
        // PASS 1 -- count every component pair's surviving points.  Pure arithmetic, nothing stored.
        auto countShellPair=[&](size_t sp)
        {
            const Shell& si=shells[sprs[sp].first]; const Shell& sj=shells[sprs[sp].second];
            // The level is a SHELL-pair property (PairLevel reads MaxExponent(i), i.e. the shared radial),
            // so compute it once and assert it per component.  Set for EVERY pair -- folded-out ones
            // included -- exactly as the per-pair build did.
            const size_t L=PairLevel(si.begin,sj.begin,ecut_L,absRelCutoff,0.0,relFieldSharp);
            for (size_t i=si.begin;i<si.end;i++)
                for (size_t j=std::max(i,sj.begin);j<sj.end;j++)
                {
                    assert(PairLevel(i,j,ecut_L,absRelCutoff,0.0,relFieldSharp)==L
                           && "the pair->level assignment must be a shell-pair property");
                    c.pairs[i*n+j].level=L;
                }
            std::vector<std::pair<size_t,size_t>> want;
            ownedPairs(sp,/*tiered*/false,want);
            if (want.empty()) return;
            const ivec3_t N=N_L[L];
            std::vector<std::pair<size_t,size_t>> here;
            ForImageOffsets(si.begin,sj.begin,A,[&](const ivec3_t& nn, const rvec3_t& Roff)
            {
                here.clear();
                for (const auto& [i,j] : want) if (!offSkip(i,j,nn)) here.push_back({i,j});
                if (here.empty()) return;                     // reduced build: member offsets never built
                ForShellPairBox(si.begin,si.end-si.begin,sj.begin,sj.end-sj.begin,Roff,A,N,
                                [&](size_t, const double* fI, const double* fJ)
                {
                    for (const auto& [i,j] : here)
                        if (std::fabs(fI[i-si.begin]*fJ[j-sj.begin])>=eps) ptsAll[i*n+j]++;
                }, eps);
            });
        };
        // PASS 2 -- materialise the streams of the pairs the budget walk tiered, and only those.
        auto buildShellPair=[&](size_t sp)
        {
            const Shell& si=shells[sprs[sp].first]; const Shell& sj=shells[sprs[sp].second];
            std::vector<std::pair<size_t,size_t>> want;
            ownedPairs(sp,/*tiered*/true,want);
            if (want.empty()) return;
            const ivec3_t N=N_L[PairLevel(si.begin,sj.begin,ecut_L,absRelCutoff,0.0,relFieldSharp)];
            std::vector<std::pair<size_t,size_t>> here;
            std::vector<PairOffsetStream> st;
            ForImageOffsets(si.begin,sj.begin,A,[&](const ivec3_t& nn, const rvec3_t& Roff)
            {
                here.clear();
                for (const auto& [i,j] : want) if (!offSkip(i,j,nn)) here.push_back({i,j});
                if (here.empty()) return;
                st.assign(here.size(),PairOffsetStream());
                for (auto& s : st) s.n=nn;
                ForShellPairBox(si.begin,si.end-si.begin,sj.begin,sj.end-sj.begin,Roff,A,N,
                                [&](size_t idx, const double* fI, const double* fJ)
                {
                    for (size_t k=0;k<here.size();k++)
                    {
                        const double val=fI[here[k].first-si.begin]*fJ[here[k].second-sj.begin];
                        if (std::fabs(val)<eps) continue;
                        st[k].Append(idx,val,tierF32[here[k].first*n+here[k].second]!=0);
                    }
                }, eps);
                for (size_t k=0;k<here.size();k++)
                    if (st[k].Points())
                        c.pairs[here[k].first*n+here[k].second].offsets.push_back(std::move(st[k]));
            });
            for (const auto& [i,j] : want) c.pairs[i*n+j].cached=true;
        };
        // One driver for both evaluation passes: shell-major, optionally threaded (the two passes have the
        // identical partition, so this is written once).
        auto runPass=[&](auto&& body)
        {
#ifdef QCHEM_OPENMP
            if (const int nthreads=PairThreads(); nthreads>1)
            {
                #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
                for (size_t sp=0;sp<sprs.size();sp++)
                {
                    try { body(sp); }
                    catch (...)
                    {
                        #pragma omp critical (gpw_pair_throw)
                        if (!firstEx) firstEx=std::current_exception();
                    }
                }
                if (firstEx) std::rethrow_exception(firstEx);
                return;
            }
#endif
            for (size_t sp=0;sp<sprs.size();sp++) body(sp);
        };

        runPass(countShellPair);
        // THE BUDGET WALK: row-major over COMPONENT pairs -- the tiering order every cache was ever built
        // in.  fp64 while the primary budget holds (bit-identical replay), fp32 for the overflow tier, drop
        // only past BOTH.  Per-pair skip-and-continue, never a lockout -- a later (smaller) pair may still
        // fit a tier.  (The old budget=0 lockout un-cached every subsequent pair after the first oversized
        // one -- the dominant cost in the 2026-07-15 multi-k profile.)  The both-tiers-exhausted "+1"
        // bookkeeping quirk is kept verbatim: the readouts report it.
        for (const auto& [i,j] : prs)
        {
            nPairs++;
            if (pairSkip(i,j)) continue;                   // fold image/dead: holds NO streams -- never
                                                           //   counted into a tier (fp64/fp32 report the
                                                           //   pairs that actually STORE)
            if (budget64==0 && budget32==0) { ptsDropped++; continue; }
            const size_t pts=ptsAll[i*n+j];
            if      (pts<=budget64) { budget64-=pts; pts64+=pts; nCached64++; tierBuild[i*n+j]=1; }
            else if (pts<=budget32) { budget32-=pts; pts32+=pts; nCached32++; tierBuild[i*n+j]=1;
                                                                              tierF32  [i*n+j]=1; }
            else                    { ptsDropped+=pts; }   // past both tiers: re-evaluated on the fly
        }
        runPass(buildShellPair);
        for (const auto& ps : c.pairs) nOffs+=ps.offsets.size();   // kept (pair, offset) terms
        c.droppedPts=ptsDropped;
        // Coverage readout per cache build (static setup, not per-iteration): the budget headroom is THE lever
        // for the analytic path's speed, so make coverage visible (dropped pairs re-evaluate every iteration).
        // Fold it into the run report's grids.stream when a run is open (it is grid metadata -- the density-
        // collocation coverage over the grid ladder); else keep the legacy one-line cerr (unbracketed builds).
        size_t nRepPairs=0;
        if (bsf) for (auto i:indices()) for (auto j:indices(i)) if (!pairSkip(i,j)) nRepPairs++;
        // The RUN count is the encoding's payoff, so report it: bytes/pt = valueWidth + 8*runs/pts, i.e.
        // the mean run length is what turns the index stream into a rounding error (see PairOffsetStream).
        size_t nRuns=0;
        for (const auto& ps : c.pairs) for (const auto& st : ps.offsets) nRuns+=st.runBase.size();
        const double meanRun = nRuns ? double(pts64+pts32)/double(nRuns) : 0.0;
        // Step 0b: announce the T3 pair-stream fold in the shared [fold] format -- ALWAYS, including when
        // it is off, because a disarmed run looks identical to an armed one on the console otherwise.
        // Armed by default on an imposed Γ run since T3.5; a NONE therefore means a FREE run, a multi-k run
        // (T3.4b) or an explicit opt-out.
        qchem::report::EmitFold("collocation streams (T3 pairs)", bsf ? bsf->maps.size() : 0,
                                nPairs, bsf ? nRepPairs : nPairs,
                                bsf ? std::string() : std::string("free/multi-k run, or GPW_STREAM_FOLD=0"));
        if (qchem::report::Depth() > 0)
        {
            qchem::report::json s;
            s["rule"]      = (absRelCutoff > 0.0 ? "abs" : "rel");
            s["relCutoff"] = absRelCutoff;
            s["pairs"]     = (long)nPairs;
            if (bsf) { s["foldOps"]=(long)bsf->maps.size(); s["repPairs"]=(long)nRepPairs; }
            s["fp64"]      = (long)nCached64;  s["pts64"] = (long)pts64;
            s["fp32"]      = (long)nCached32;  s["pts32"] = (long)pts32;
            s["offsets"]   = (long)nOffs;      // kept (pair, image-offset) terms across all pairs
            s["runs"]      = (long)nRuns;      // contiguous index runs (the stored geometry)
            s["meanRun"]   = meanRun;          //   and their mean length -- the bandwidth lever
            s["dropped"]   = (long)ptsDropped;
            s["budget64"]  = (long)BudgetPts(); s["budget32"] = (long)BudgetPtsF32();
            qchem::report::EmitAt("grids", "stream", s);
        }
        else
        {
            std::cerr << "[stream cache] shape=(";
            for (size_t l=0;l<N_L.size();l++) std::cerr << (l?",":"") << N_L[l].x;
            std::cerr << ") rule=" << (absRelCutoff>0.0 ? "abs " : "rel ") << absRelCutoff;
            if (bsf) std::cerr << "  FOLD |ops|=" << bsf->maps.size() << " repPairs " << nRepPairs << "/" << nPairs;
            std::cerr << "  pairs " << nPairs
                      << ": fp64 " << nCached64 << " (" << pts64 << " pts), fp32 " << nCached32
                      << " (" << pts32 << " pts), offsets " << nOffs
                      << ", runs " << nRuns << " (mean " << meanRun << " pts), dropped " << ptsDropped
                      << " pts (budgets " << BudgetPts() << "/" << BudgetPtsF32() << ")" << std::endl;
        }
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
                                         const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L,
                                         double relFieldSharp=-1.0) const
    {
        const size_t K=N_L.size();
        assert(K>0 && ecut_L.size()==K);
        // THE per-iteration cost bucket (doc/GPWPlan1.md "fast-recompute campaign"): every SCF step
        // scatters every (pair, offset) term, either as a cached-stream replay or -- past the budget --
        // as a live box evaluation.  The stream BUILD nested inside (EnsureStreams) charges its own
        // bucket, so this one reads as the pure per-iteration scatter.
        qchem::report::Timed timer("scf: collocate density (pair scatter)");
        std::vector<rvec_t> rho(K);
        for (size_t l=0; l<K; l++) rho[l]=rvec_t(size_t(N_L[l].x)*N_L[l].y*N_L[l].z, 0.0);
        // Hermitian fold: the (j,i,-R) term contributes the SAME (idx, value, weight) as (i,j,R) -- the product
        // field is the lattice-translated twin and D_ji e^{-ik(-R)} = conj(D_ij e^{-ikR}) -- so loop j>=i and
        // double the off-diagonal weight (a 2x saving over the full n^2 loop, exact).
        const StreamCache& sc=EnsureStreams(A,N_L,ecut_L,0.0,relFieldSharp);   // density: the relative smooth-field rule
        const size_t nn=size();
        const StreamFold* sf=itsStreamFold.get();               // T3 route (b): reduced scatter (§6b), null = full
        // T3.5: the reduced scatter reads the ORBIT-PROJECTED D at each representative, never the
        // representative's own element -- see FoldProjectedD.  O(n^2) per call beside an O(pairs x points)
        // scatter, and it is what makes the folded density equal the unfolded imposed one for ANY iterate.
        const std::vector<dcmplx> Dsym=FoldProjectedD(D);
        // TWO PASSES over the pair set (stage B, 2026-08-19).  A pair with a CACHED stream replays it here,
        // in the original ROW-MAJOR order, so a fully-cached run (every anchor) accumulates into `rho` in
        // exactly the order it always did -- bit-identical.  Pairs with NO stream (over the budget) are left
        // to `scatterShell` below, which walks their boxes SHELL-BLOCKED.  The split is what keeps the
        // in-budget regime untouched while the over-budget one gets the kernel.
        //
        // Per-pair scatter into a caller-supplied set of level densities (`dst`): the SAME `rho` when serial,
        // a private per-thread accumulator when threaded -- so the arithmetic PER PAIR is bit-identical either
        // way; only the cross-pair reduction order changes (see PairThreads).
        auto scatter=[&](size_t i, size_t j, std::vector<rvec_t>& dst)
        {
            // Reduced mode: only orbit-REPRESENTATIVE pairs scatter, each rep offset weighted by its full
            // triple-orbit multiplicity (pairMult x within) and the ORBIT-PROJECTED D (T3.5); the caller's
            // dense group-average completes it.  The D-aware kill stays on the MEMBER weight |c|
            // (orbit-invariant -- exactly so now that the projected D is one value per orbit), so reduced
            // and full keep the identical active set (§6b item 6).
            const std::map<long long,unsigned>* fm=nullptr;
            double pmul=1.0;
            dcmplx Dij(D(i,j));
            if (sf)
            {
                const PairEdge& pe=sf->pairs[i*nn+j];
                if (pe.dead || pe.rep!=int(i*nn+j)) return;     // image/dead pair: its rep carries the orbit
                fm=&sf->offMult[i*nn+j];
                pmul=double(pe.pairMult);
                Dij=Dsym[i*nn+j];
            }
            if (Dij==dcmplx(0.0)) return;
            const double fold=(i==j)?1.0:2.0;
            const PairStreams& ps=sc.pairs[i*nn+j];
            rvec_t& r=dst[ps.level];
            if (ps.cached)
                for (const PairOffsetStream& st : ps.offsets)
                {
                    double wm=1.0;
                    if (fm)
                    {
                        const auto it=fm->find(OffKey(st.n));
                        assert(it!=fm->end());
                        if (it==fm->end() || it->second==0) continue;   // member offset: its rep carries it
                        wm=pmul*double(it->second);
                    }
                    const double c=fold*std::real(Dij*std::conj(phase(st.n)));  // Re[D_ij e^{-ik.R_n}] offset weight
                    if (std::fabs(c)*st.maxv < kDensityEps()) continue;   // D-aware kill: sub-eps density term
                    const double cw=c*wm;
                    // Run-major replay: one CONTIGUOUS destination span per run, values streamed in build
                    // order -- the same points in the same order as the per-point-index form, so the sums
                    // are bit-identical; the index stream is gone (see PairOffsetStream).
                    const unsigned* rb=st.runBase.data();
                    const unsigned* rl=st.runLen.data();
                    const size_t nr=st.runBase.size();
                    if (!st.val.empty())
                    {
                        const double* v=st.val.data();
                        for (size_t q=0, k=0; q<nr; q++)
                        {
                            double* d=&r[rb[q]];
                            for (unsigned t=0, m=rl[q]; t<m; t++, k++) d[t]+=cw*v[k];
                        }
                    }
                    else                                            // fp32 overflow tier (values demoted)
                    {
                        const float* v=st.val32.data();
                        for (size_t q=0, k=0; q<nr; q++)
                        {
                            double* d=&r[rb[q]];
                            for (unsigned t=0, m=rl[q]; t<m; t++, k++) d[t]+=cw*double(v[k]);
                        }
                    }
                }
            // (no `else`: an uncached pair is scattered SHELL-BLOCKED by scatterShell below)
        };
        // The SHELL-BLOCKED on-the-fly scatter (stage B): the over-budget pairs of ONE shell pair, whose box
        // walk -- offsets, level, geometry, ellipsoid pre-screen, wrap -- is shared by all of them.
        //
        // THE D-AWARE BOX, resolved across the shell pair (the design decision, doc/GPWPlan1.md "Round 4"):
        // per component pair the tolerance is still eps/|c_ij| (its own |D| weight), but ONE box must serve
        // the whole shell pair, so the box is sized from the UNION -- the TIGHTEST eps present, i.e. the
        // LARGEST box -- and each component pair then applies its OWN |val| < eps_ij magnitude screen.  That
        // only ever ADDS sub-eps terms and never drops one, so it keeps the no-cut discipline; the price is
        // that one large-|D| component makes its siblings pay the bigger box.  (Terms a per-component box
        // would have killed WHOLE via the prefactor early-out are instead killed point-by-point -- same
        // result, a little more work.)
        const std::vector<Shell>& shells=Shells();
        std::vector<std::pair<size_t,size_t>> sprs;
        for (size_t a=0;a<shells.size();a++) for (size_t b=a;b<shells.size();b++) sprs.push_back({a,b});
        auto scatterShell=[&](size_t sp, std::vector<rvec_t>& dst)
        {
            const Shell& si=shells[sprs[sp].first]; const Shell& sj=shells[sprs[sp].second];
            // The component pairs of this shell pair with NO stream, a live D and (under a fold) a
            // representative edge -- exactly the pairs `scatter` above declined.
            std::vector<size_t> want;                                   // encoded i*nn+j
            for (size_t i=si.begin;i<si.end;i++)
                for (size_t j=std::max(i,sj.begin);j<sj.end;j++)
                {
                    if (sc.pairs[i*nn+j].cached) continue;
                    if (sf) { const PairEdge& pe=sf->pairs[i*nn+j];
                              if (pe.dead || pe.rep!=int(i*nn+j)) continue; }
                    if ((sf ? Dsym[i*nn+j] : dcmplx(D(i,j)))==dcmplx(0.0)) continue;
                    want.push_back(i*nn+j);
                }
            if (want.empty()) return;
            const size_t L=sc.pairs[want.front()].level;                // a shell-pair property (see the build)
            rvec_t& r=dst[L];
            std::vector<size_t> here;                                   // this offset's live pairs
            std::vector<double> cwHere, epsHere;
            ForImageOffsets(si.begin,sj.begin,A,[&](const ivec3_t& n, const rvec3_t& Roff)
            {
                here.clear(); cwHere.clear(); epsHere.clear();
                double epsUnion=0.0;
                const double pf=PairPrefactorExp(si.begin,sj.begin,Roff);   // screen (1), shared by the shell pair
                for (size_t key : want)
                {
                    const size_t i=key/nn, j=key%nn;
                    double wm=1.0;
                    if (sf)
                    {
                        const auto& mm=sf->offMult[key];
                        const auto it=mm.find(OffKey(n));
                        assert(it!=mm.end());
                        if (it==mm.end() || it->second==0) continue;    // member offset: its rep carries it
                        wm=double(sf->pairs[key].pairMult)*double(it->second);
                    }
                    const dcmplx Dij = sf ? Dsym[key] : dcmplx(D(i,j));        // T3.5: the orbit-projected D
                    const double c=((i==j)?1.0:2.0)*std::real(Dij*std::conj(phase(n)));
                    if (c==0.0) continue;
                    // MEMBER rule (|c|, not |c*wm|) so the reduced box matches every orbit member's.
                    const double e=std::max(kDensityEps(), kDensityEps()/std::fabs(c));
                    if (pf>=-std::log(e)) continue;                     // whole box below THIS pair's eps
                    here.push_back(key); cwHere.push_back(c*wm); epsHere.push_back(e);
                    epsUnion = epsUnion==0.0 ? e : std::min(epsUnion,e);
                }
                if (here.empty()) return;
                ForShellPairBox(si.begin,si.end-si.begin,sj.begin,sj.end-sj.begin,Roff,A,N_L[L],
                                [&](size_t idx, const double* fI, const double* fJ)
                {
                    for (size_t k=0;k<here.size();k++)
                    {
                        const double v=fI[here[k]/nn-si.begin]*fJ[here[k]%nn-sj.begin];
                        if (std::fabs(v)<epsHere[k]) continue;          // each pair keeps its OWN screen
                        r[idx]+=cwHere[k]*v;
                    }
                }, epsUnion);
            });
        };
#ifdef QCHEM_OPENMP
        if (const int nthreads=PairThreads(); nthreads>1)
        {
            // OpenMP over the pairs.  Each thread owns a PRIVATE level-density accumulator `mine` (declared
            // inside the parallel region -> genuinely thread-local, so no omp_get_thread_num / <omp.h>
            // needed), scatters its share of the pairs into it, then folds it into the shared rho under a
            // critical section (the "per-thread accumulators + reduce").  The reduction reorders the
            // cross-pair grid sums, so a threaded run drifts a few ULPs from serial -- accepted; the Si
            // bit-anchors always run serial (GPW_OMP_THREADS unset -> nthreads==1 -> the branch below).
            std::vector<std::pair<size_t,size_t>> prs;
            for (auto i:indices()) for (auto j:indices(i)) prs.push_back({i,j});
            // THROW CONTAINMENT: an exception escaping an OpenMP construct is an instant std::terminate
            // (the intermittent threaded-NaF abort: "terminate called recursively" -- several workers
            // throwing at once).  Capture the FIRST worker exception, finish the region, rethrow it
            // serially -- so a throw surfaces as the ordinary exception it is on the serial path.
            std::exception_ptr firstEx;
            #pragma omp parallel num_threads(nthreads)
            {
                std::vector<rvec_t> mine(K);
                for (size_t l=0;l<K;l++) mine[l]=rvec_t(rho[l].size(),0.0);
                #pragma omp for schedule(dynamic) nowait
                for (size_t p=0;p<prs.size();p++)
                {
                    try { scatter(prs[p].first,prs[p].second,mine); }
                    catch (...)
                    {
                        #pragma omp critical (gpw_pair_throw)
                        if (!firstEx) firstEx=std::current_exception();
                    }
                }
                #pragma omp for schedule(dynamic) nowait
                for (size_t sp=0;sp<sprs.size();sp++)          // the over-budget remainder, shell-blocked
                {
                    try { scatterShell(sp,mine); }
                    catch (...)
                    {
                        #pragma omp critical (gpw_pair_throw)
                        if (!firstEx) firstEx=std::current_exception();
                    }
                }
                #pragma omp critical
                for (size_t l=0;l<K;l++)
                    for (size_t g=0, m=rho[l].size(); g<m; g++) rho[l][g]+=mine[l][g];
            }
            if (firstEx) std::rethrow_exception(firstEx);
            return rho;
        }
#endif
        for (auto i:indices()) for (auto j:indices(i)) scatter(i,j,rho);   // serial (byte-identical default)
        for (size_t sp=0;sp<sprs.size();sp++) scatterShell(sp,rho);       // the over-budget remainder
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
        double absRelCutoff=0.0;                                //   + the assignment rule key (0=relative, >0=absolute Ha/exponent)
        double relFieldSharp=-1.0;                              //   + the relative rule's beta floor (level-assignment identity)
        std::vector<rvec_t>  V_L;                               //   + the integrated field itself (exact)
        std::vector<PairB>   B;                                 // indexed [i*n+j], j>=i only
    };
    static constexpr size_t kMaxIntegrateMemos=4;               // static PP + the iteration's KS fields
    bool SameShape(const IntegrateMemo& m, const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L,
                   double absRelCutoff, double relFieldSharp) const
    {
        if (m.absRelCutoff!=absRelCutoff || m.relFieldSharp!=relFieldSharp
         || m.ecut_L!=ecut_L || m.N_L.size()!=N_L.size()) return false;
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
                               double absRelCutoff=0.0, const chmat_t* screenD=nullptr,
                               double fieldSharpness=0.0,         // beta_loc: the sharp field's own exponent (local-PP)
                               double relFieldSharp=-1.0,         // the relative rule's beta floor (must match collocation)
                               const std::vector<size_t>* pairLevels=nullptr) const  // explicit assignment (static fields)
    {
        const size_t K=N_L.size();
        assert(K>0 && ecut_L.size()==K && V_L.size()==K);
        assert(!pairLevels || pairLevels->size()==size()*size());
        const size_t nn=size();
        // The collocation's sibling per-iteration bucket -- charged ONLY on the density-rule path (the KS
        // field, once per k-block per iteration).  The STATIC sharp-field sweeps (local PP, explicit
        // pairLevels) stay uncharged on purpose: they already sit inside their own named setup buckets
        // (setup: local-PP LONG/SHORT, separable PP), and an exclusive child here would hollow those out.
        std::optional<qchem::report::Timed> timer;
        if (absRelCutoff==0.0 && !pairLevels) timer.emplace("scf: integrate-back (pair gather)");
        chmat_t h(size());
        // T3 route (b) (§6b): the reduced gather runs only on the per-iteration density-rule path (the
        // static sharp-field sweeps and explicit-level calls keep their own machinery, incl. T1).  It
        // bypasses the memo like the screened calls: the reduced sweep is cheap by construction.
        const StreamFold* sf = (absRelCutoff==0.0 && !pairLevels) ? itsStreamFold.get() : nullptr;
        const bool memoize = (screenD==nullptr && !sf);
        // T3.5: under a fold the D-aware screen must be ORBIT-INVARIANT (every member of a pair orbit
        // truncated the same way, since one representative now stands for all of them) -- the orbit MAX of
        // the screen magnitude.  See FoldScreenMax for why this is NOT the signed FoldProjectedD.
        const std::vector<double> Dscreen = (sf && screenD) ? FoldScreenMax(*screenD) : std::vector<double>{};
        // Memo hit: the per-offset reductions are already computed for this exact field -- contract phases only.
        if (memoize)
            for (const auto& m : itsIntegrateMemos)
                if (SameShape(m,N_L,ecut_L,absRelCutoff,relFieldSharp) && SameField(m.V_L,V_L))
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
        // absRelCutoff > 0 marks a STATIC sharp-field call (the local PP, built once per SCF): evaluate on
        // the fly -- caching a single-shot sweep wastes the whole budget the per-iteration path needs.  (Its
        // B reductions land in the memo, so the other k-blocks never repeat the sweep.)  An EXPLICIT
        // pairLevels assignment likewise bypasses the streams (they were built at the internal rule's
        // levels); its memo key is honest because the caller's ladder (the custom top level) already
        // distinguishes it from every internal-rule shape.
        const StreamCache* sc = (absRelCutoff==0.0 && !pairLevels)
                              ? &EnsureStreams(A,N_L,ecut_L,absRelCutoff,relFieldSharp) : nullptr;
        IntegrateMemo* memo=nullptr;
        if (memoize)
        {
            if (itsIntegrateMemos.size()>=kMaxIntegrateMemos) itsIntegrateMemos.erase(itsIntegrateMemos.begin());
            itsIntegrateMemos.emplace_back();
            memo=&itsIntegrateMemos.back();
            memo->N_L=N_L; memo->ecut_L=ecut_L; memo->absRelCutoff=absRelCutoff;
            memo->relFieldSharp=relFieldSharp; memo->V_L=V_L;
            memo->B.resize(nn*nn);
        }
        // Per-pair integrate-back.  Each pair writes ONLY its own h(i,j) and its own (pre-sized, disjoint)
        // memo->B[i*nn+j] slot -- so this is embarrassingly parallel with NO reduction (unlike collocation,
        // whose pairs scatter into a shared grid).
        auto integratePair=[&](size_t i, size_t j)           // j>=i (Hermitian upper triangle)
        {
            // Reduced mode (§6b item 5): gather REPRESENTATIVE pairs only, rep offsets weighted by their
            // within-pair multiplicity; image pairs are filled by the representation transform after the
            // loop (h_{i'j'} = sigma h_ij; dead pairs h = 0 -- the projected value).
            const std::map<long long,unsigned>* fm=nullptr;
            if (sf)
            {
                const PairEdge& pe=sf->pairs[i*nn+j];
                if (pe.dead || pe.rep!=int(i*nn+j)) return;
                fm=&sf->offMult[i*nn+j];
            }
            const size_t  l = pairLevels ? (*pairLevels)[i*nn+j]
                            : sc         ? sc->pairs[i*nn+j].level
                            :              PairLevel(i,j,ecut_L,absRelCutoff,fieldSharpness,relFieldSharp);
            const rvec_t& V=V_L[l];
            const double  w=A.GetCellVolume()/double(V.size());   // the level's quadrature weight Omega/Npts(l)
            // pb.nb records the per-offset B(n) reductions FOR THE MEMO only (the phase-independent replay).  On
            // the per-iteration density path (memo==null, V changes each SCF step so the memo can never hit) the
            // recording is pure waste -- pb points at a throwaway `dummy`, and emplace_back'ing ~1M offsets into
            // it per run just churns allocations.  Record only when there's a memo to serve.
            PairB  dummy;
            PairB& pb = memo ? memo->B[i*nn+j] : dummy;
            if (memo) pb.level=l;
            // The D-aware offset screen.  It used to reuse the density collocation's COEFFICIENT,
            // fold*Re[D_ij e^{-ik.R_n}], so that both directions kept the identical active set.  That is
            // correct on the COLLOCATION side -- there Re[D e^{-ikR}] multiplies a real pair product, so a
            // zero coefficient contributes exactly nothing to rho -- and WRONG here, because a real part is
            // not a magnitude.  This direction's term is phase(n)*b, whose size is |b| however the phase is
            // oriented, so a term with Re[D conj(phase)] = 0 still adds a purely IMAGINARY, entirely
            // non-negligible amount to H.
            //
            // At a QUARTER-INTEGER k that is not an accident but a systematic kill: e^{2 pi i k n} = i^n is
            // purely imaginary for every ODD offset, so for real-ish D the old test discarded EVERY odd-offset
            // term and the Hartree/XC matrix came out exactly real (measured: maxIm(dV) = 0 at k=1/4 against
            // 0.067 at k=0.25001).  The resulting H is missing its imaginary part, its spectrum is wrong, and
            // the SCF converges to a state 2.5 Ha high -- doc/Benchmark.md footnote 1, the shifted-MP defect.
            // The screen is now the true magnitude |D_ij| (= |D_ij conj(phase)|, since |phase|=1), which is
            // what "the magnitude screen is the only truncation" means; it is strictly MORE conservative than
            // the old test, so it can only keep terms the old one dropped.
            const double fold=(i==j)?1.0:2.0;
            const double Dmag = !screenD          ? 0.0
                              : Dscreen.empty()   ? fold*std::abs(dcmplx((*screenD)(i,j)))
                              :                     fold*Dscreen[i*nn+j];
            dcmplx s(0.0);
            if (sc && sc->pairs[i*nn+j].cached)
                for (const PairOffsetStream& st : sc->pairs[i*nn+j].offsets)
                {
                    double wm=1.0;
                    if (fm)
                    {
                        const auto it=fm->find(OffKey(st.n));
                        assert(it!=fm->end());
                        if (it==fm->end() || it->second==0) continue;   // member offset: its rep carries it
                        wm=double(it->second);
                    }
                    if (screenD && Dmag*st.maxv < kDensityEps()) continue;   // TRUE magnitude, not Re[]
                    double b=0.0;
                    // Run-major gather -- the exact mirror of the collocation scatter (same runs, same
                    // order), so the reduction stays bit-identical to the per-point-index form and the
                    // collocate/integrate adjoint stays machine-exact.
                    const unsigned* rb=st.runBase.data();
                    const unsigned* rl=st.runLen.data();
                    const size_t nr=st.runBase.size();
                    if (!st.val.empty())
                    {
                        const double* v=st.val.data();
                        for (size_t q=0, k=0; q<nr; q++)
                        {
                            const double* Vp=&V[rb[q]];
                            for (unsigned t=0, m=rl[q]; t<m; t++, k++) b+=v[k]*Vp[t];
                        }
                    }
                    else                                            // fp32 overflow tier (values demoted)
                    {
                        const float* v=st.val32.data();
                        for (size_t q=0, k=0; q<nr; q++)
                        {
                            const double* Vp=&V[rb[q]];
                            for (unsigned t=0, m=rl[q]; t<m; t++, k++) b+=double(v[k])*Vp[t];
                        }
                    }
                    if (memo) pb.nb.emplace_back(st.n,b);
                    s+=phase(st.n)*(wm*b);
                }
            else return;      // static sharp-field call or over-budget pair: gathered SHELL-BLOCKED below
            s*=w;
            h(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) auto-set to conj
        };
        // The SHELL-BLOCKED on-the-fly gather (stage B) -- the exact mirror of `scatterShell`.  It serves the
        // over-budget pairs AND, since those calls hold no stream cache at all, the whole STATIC SHARP-FIELD
        // sweep (MakeLocalPP's kappa rule, the explicit-pairLevels V_loc ball): every pair of a shell pair
        // shares one box walk.  The D-aware box resolves as the UNION with per-component screens, exactly as
        // on the density side, so collocate and integrate keep the SAME active set and stay adjoints.
        // Each pair still writes only its own h(i,j) and its own memo slot -- no reduction, so unlike the
        // scatter this pass may run in any order.
        const std::vector<Shell>& shells=Shells();
        std::vector<std::pair<size_t,size_t>> sprs;
        for (size_t a=0;a<shells.size();a++) for (size_t b=a;b<shells.size();b++) sprs.push_back({a,b});
        auto integrateShell=[&](size_t sp)
        {
            const Shell& si=shells[sprs[sp].first]; const Shell& sj=shells[sprs[sp].second];
            std::vector<size_t> want;                                   // encoded i*nn+j
            for (size_t i=si.begin;i<si.end;i++)
                for (size_t j=std::max(i,sj.begin);j<sj.end;j++)
                {
                    if (sc && sc->pairs[i*nn+j].cached) continue;       // replayed by integratePair
                    if (sf) { const PairEdge& pe=sf->pairs[i*nn+j];
                              if (pe.dead || pe.rep!=int(i*nn+j)) continue; }
                    want.push_back(i*nn+j);
                }
            if (want.empty()) return;
            // The level is a shell-pair property in all three regimes: the explicit pairLevels assignment
            // (StaticFieldPairLevels keys on MaxExponent(i)+MaxExponent(j)), the stream cache's stored
            // level, and PairLevel itself all read the shared radial.
            const size_t i0=want.front()/nn, j0=want.front()%nn;
            const size_t l = pairLevels ? (*pairLevels)[i0*nn+j0]
                           : sc         ? sc->pairs[i0*nn+j0].level
                           :              PairLevel(i0,j0,ecut_L,absRelCutoff,fieldSharpness,relFieldSharp);
            const rvec_t& V=V_L[l];
            const double  w=A.GetCellVolume()/double(V.size());
#ifndef NDEBUG
            for (size_t key : want)                              // every component pair must agree on it
                assert(l==(pairLevels ? (*pairLevels)[key]
                          : sc        ? sc->pairs[key].level
                          :             PairLevel(key/nn,key%nn,ecut_L,absRelCutoff,fieldSharpness,relFieldSharp))
                       && "the pair->level assignment must be a shell-pair property");
#endif
            std::vector<dcmplx> s(want.size(), dcmplx(0.0));      // (integratePair already set memo->B[].level)
            std::vector<size_t> here;
            std::vector<double> wmHere, epsHere, bHere;
            ForImageOffsets(si.begin,sj.begin,A,[&](const ivec3_t& n, const rvec3_t& Roff)
            {
                here.clear(); wmHere.clear(); epsHere.clear();
                double epsUnion=0.0;
                const double pf=PairPrefactorExp(si.begin,sj.begin,Roff);   // screen (1), shared by the shell pair
                for (size_t k=0;k<want.size();k++)
                {
                    const size_t key=want[k], i=key/nn, j=key%nn;
                    double wm=1.0;
                    if (sf)
                    {
                        const auto& mm=sf->offMult[key];
                        const auto it=mm.find(OffKey(n));
                        assert(it!=mm.end());
                        if (it==mm.end() || it->second==0) continue;
                        wm=double(it->second);
                    }
                    double e=kDensityEps();      // collocation floor (decoupled from the analytic kScreenEps)
                    if (screenD)
                    {
                        const double c=std::fabs(((i==j)?1.0:2.0)
                                                 *std::real(dcmplx((*screenD)(i,j))*std::conj(phase(n))));
                        if (c==0.0) continue;
                        e=std::max(kDensityEps(), kDensityEps()/c);
                    }
                    if (pf>=-std::log(e)) continue;                     // whole box below THIS pair's eps
                    here.push_back(k); wmHere.push_back(wm); epsHere.push_back(e);
                    epsUnion = epsUnion==0.0 ? e : std::min(epsUnion,e);
                }
                if (here.empty()) return;
                bHere.assign(here.size(),0.0);
                ForShellPairBox(si.begin,si.end-si.begin,sj.begin,sj.end-sj.begin,Roff,A,N_L[l],
                                [&](size_t idx, const double* fI, const double* fJ)
                {
                    const double Vv=V[idx];
                    for (size_t q=0;q<here.size();q++)
                    {
                        const size_t key=want[here[q]];
                        const double v=fI[key/nn-si.begin]*fJ[key%nn-sj.begin];
                        if (std::fabs(v)<epsHere[q]) continue;          // each pair keeps its OWN screen
                        bHere[q]+=v*Vv;
                    }
                }, epsUnion);
                for (size_t q=0;q<here.size();q++)
                {
                    if (memo) memo->B[want[here[q]]].nb.emplace_back(n,bHere[q]);
                    s[here[q]]+=phase(n)*(wmHere[q]*bHere[q]);
                }
            });
            for (size_t k=0;k<want.size();k++)
            {
                const size_t i=want[k]/nn, j=want[k]%nn;
                const dcmplx sk=s[k]*w;
                h(i,j) = (i==j) ? dcmplx(std::real(sk),0.0) : sk;
            }
        };
        // Reduced mode: fill the image pairs from their reps by the representation transform (Hermitian
        // flip = conjugate; at Γ everything is real).  Dead pairs get the projected value 0.
        auto fillImages=[&]()
        {
            if (!sf) return;
            for (auto i:indices()) for (auto j:indices(i))
            {
                const PairEdge& pe=sf->pairs[i*nn+j];
                if (pe.dead) { h(i,j)=dcmplx(0.0); continue; }
                if (pe.rep==int(i*nn+j)) continue;
                const size_t a=size_t(pe.rep)/nn, b=size_t(pe.rep)%nn;
                // h_{i'j'} = zeta h_{ab} (§6b/T3.4); a flipped edge stores the Hermitian twin slot, so it
                // takes the conjugate of the WHOLE product (at Γ zeta is ±1 and conj is the identity).
                dcmplx v = pe.zeta*dcmplx(h(a,b));
                if (pe.flip) v=std::conj(v);
                h(i,j) = (i==j) ? dcmplx(std::real(v),0.0) : v;
            }
        };
#ifdef QCHEM_OPENMP
        if (const int nthreads=PairThreads(); nthreads>1)
        {
            std::vector<std::pair<size_t,size_t>> prs;
            for (auto i:indices()) for (auto j:indices(i)) prs.push_back({i,j});
            std::exception_ptr firstEx;   // throw containment -- see CollocateDensity
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (size_t p=0;p<prs.size();p++)
            {
                try { integratePair(prs[p].first,prs[p].second); }
                catch (...)
                {
                    #pragma omp critical (gpw_pair_throw)
                    if (!firstEx) firstEx=std::current_exception();
                }
            }
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (size_t sp=0;sp<sprs.size();sp++)              // the uncached remainder, shell-blocked
            {
                try { integrateShell(sp); }
                catch (...)
                {
                    #pragma omp critical (gpw_pair_throw)
                    if (!firstEx) firstEx=std::current_exception();
                }
            }
            if (firstEx) std::rethrow_exception(firstEx);
            fillImages();
            return h;
        }
#endif
        for (auto i:indices()) for (auto j:indices(i)) integratePair(i,j);   // serial (byte-identical default)
        for (size_t sp=0;sp<sprs.size();sp++) integrateShell(sp);           // the uncached remainder
        fillImages();
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
    // One monomial term of a per-atom local Gaussian operator: the atom's GaussianRF (SHARED across its
    // polynomial terms -- built once per atom image), the term's Cartesian-monomial Polarization, and its coeff.
    struct OpTerm { std::shared_ptr<GaussianRF> gr; Polarization pol; double c; };
    // Flatten the structure's per-atom operators (opForZ(Z_a) at R_a) into (gr, pol, c) terms for the 3-centre
    // Overlap3C loop.  For a PERIODIC cell (\a A non-null) the local potential is the FIXED periodic field, so
    // each atom is placed at every IMAGE within (operator reach + max orbital reach) of the home atom: an
    // operator co-located with an imaged chi_j atom that lands near a home orbital contributes, so home-cell
    // atoms alone under-count (measured: Si Gamma -7.67 vs -7.11).  The operator is compact (Gaussian tail), so
    // the image shell is finite; a far one still lands in the list but Overlap3C screens it to ~0.  \a A null =
    // the FINITE (molecule/box) limit: home atoms only.  One GaussianRF per atom image, shared by its terms.
    std::vector<OpTerm> BuildOpTerms(const Structure* cl,
                                     const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ,
                                     const UnitCell* A) const
    {
        double maxOrbReach=0.0;                              // the product of two orbitals can spread this far
        for (auto i:indices())
        {
            double amin=radials[i]->GetExponents()[0];
            for (double e:radials[i]->GetExponents()) amin=std::min(amin,e);
            maxOrbReach=std::max(maxOrbReach, std::sqrt(-std::log(kScreenEps())/amin));
        }
        std::vector<OpTerm> ops;
        for (Atom* a : *cl)
        {
            Molecule::LatticeSum1E::GaussianFunction g=opForZ(a->itsZ);
            const int Lg=GaussDegree(g);
            auto place=[&](const rvec3_t& center)
            {
                auto gr=std::make_shared<GaussianRF>(g.alpha, center, Lg);
                for (const auto& t : g.terms) ops.push_back({gr, Polarization(t.p), t.c});
            };
            if (!A) { place(a->itsR); continue; }           // finite: home atom only
            const double rEnum=std::sqrt(-std::log(kScreenEps())/g.alpha)+maxOrbReach;   // operator reach + orbital reach
            for (const auto& n : A->CellsInSphere(rEnum + A->GetMaximumCellEdge()))
            {
                const rvec3_t Roff=A->ToCartesian(rvec3_t(double(n.x),double(n.y),double(n.z)));
                if (Roff.x*Roff.x+Roff.y*Roff.y+Roff.z*Roff.z <= rEnum*rEnum) place(a->itsR+Roff);
            }
        }
        return ops;
    }

    mutable std::vector<Shell>         itsShells;          //!< the shell partition (lazy; geometry-fixed)
    mutable std::vector<StreamCache>   itsStreamCaches;    //!< analytic pair-box streams, one per ladder shape
    mutable std::vector<IntegrateMemo> itsIntegrateMemos;  //!< phase-independent B_ij(n) per exact field (LRU)
    mutable std::unique_ptr<StreamFold> itsStreamFold;     //!< T3 route (b) (pair, R) orbit fold (null = free run)
};

static_assert(is1E_Evaluator <NR_Evaluator>, "NR_Evaluator must satisfy is1E_Evaluator");
static_assert(isDFT_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isDFT_Evaluator (ThreeC kernel)");
static_assert( isHF_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isHF_Evaluator (FourC kernel)");

} //namespace

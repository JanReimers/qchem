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
#include <iostream>   // std::cout (the one-line [collocation] kernel + task-list provenance)
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
#include <cstring>    // std::memcpy (the integrate-census field hash)
// GPW pair-loop parallelism (opt-in via GPW_OMP_THREADS).  Gated on QCHEM_OPENMP -- our OWN macro (defined
// by CMake when the OpenMP flags are applied) rather than _OPENMP, because this LLVM toolchain ships no
// libomp and the libgomp fallback (-fopenmp=libgomp) honours the pragmas but does NOT define _OPENMP.  The
// private-buffer + critical-reduce pattern below needs no <omp.h> (no omp_*() calls), only the pragmas.
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;      // PGData
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;  // GaussianRF named kernels
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;// Polarization
import qchem.IntPow;                                             // uintpow (the monomial power tables of the collocation box walk)
import qchem.BasisSet.Molecule.LatticeSum1E;                       // GaussianFunction (the <chi_i|g> seam type)
export import qchem.Symmetry.Lattice_3D.SpaceGroup;                // DirectOp {W|τ} (the T3 stream-fold ops)
import qchem.Structure;
import qchem.UnitCell;                                             // UnitCell (ToCartesian/ToFractional: grid<->cell for collocation)
import qchem.Types;
import qchem.Blaze;                                                // rsmat_t (the lattice-sum matrices)
import qchem.Reporting;                                            // fold the task-list geometry into grids.boxTasks

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
    //! \brief \c GPW_DAWARE_SCREEN=0 makes the box tolerance GEOMETRY-ONLY: \f$\varepsilon\f$ flat, no
    //! \f$\varepsilon/|c_{ij}|\f$ widening and no per-(pair, offset) \f$D\f$ kill.  ⚠ AN EXPERIMENT'S
    //! KNOB, NOT THE DESIGN.  The design this is measuring the case for is a \c LatticeScreener interface
    //! with \c DAwareScreener / \c GeometryOnlyScreener behind it, chosen once at the \c SolidCalculation
    //! level — see doc/ScreeningPlan.md.  Do not grow more branches on this bool; replace it.
    //!
    //! WHY IT IS WORTH MEASURING (2026-08-28): CP2K does NOT screen on the density
    //! (\c task_list_methods.F: the radius takes no density argument, its eps is the global
    //! \c eps_rho_rspace, and \c radius_list is fixed task-list data), and
    //! \f$\varepsilon_{eff}\f$ is the ONLY \f$D\f$-dependent input to \c BoxGeom — so switching it off
    //! makes the whole box geometry ITERATION-INVARIANT and hoistable into the task list.
    static bool UseDAwareScreen()
    { static const bool b=[]{const char* s=std::getenv("GPW_DAWARE_SCREEN"); return !s || std::atoi(s)!=0;}(); return b; }

    //! The lattice-sum ECONOMY report (LatticeSum1E face; doc there).  Reach numbers use the SAME formulas
    //! as the screens: a single orbital's tail reach is sqrt(-ln eps/alpha), a pair's conservative reach is
    //! the sum -- worst pair = diffuse x diffuse = 2*sqrt(-ln eps/alpha_min).  The CellsInSphere count over
    //! that reach is the per-pair image-enumeration size a grad student sees JUMP when diffuse functions
    //! arrive; the screens then prune it per term (kept counts = the [collocation] task-list readout).
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
                 <<" cells (kept counts on the [collocation] line)"<<std::endl;
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
    //! Fill \a pw[0..L] with \f$x^k\f$ using EXACTLY \c uintpow's associations -- including its \f$x=0\f$
    //! and \f$k=0\f$ conventions -- so indexing this table is BIT-IDENTICAL to the \c uintpow call it
    //! replaces.  \a L is bounded by \c kMaxDeg (the caller sizes \a pw).
    static void PowerTable(double x, int L, double* pw)
    {
        pw[0]=1.0;
        if (L<1) return;
        if (x==0.0) { for (int k=1;k<=L;k++) pw[k]=0.0; return; }   // uintpow's x==0 early-out (+0.0, never -0.0)
        pw[1]=x;
        if (L>=2) pw[2]=x*x;
        if (L>=3) pw[3]=x*x*x;
        if (L>=4) pw[4]=(x*x)*(x*x);
        for (int k=5;k<=L;k++) pw[k]=uintpow(x,(unsigned)k);        // the binary-powering tail, unchanged
    }

    //! \brief The BOX GEOMETRY of one (shell pair, cross-cell offset): everything the walk derives before
    //! it touches a grid point.  Extracted 2026-08-27 so the walk and the separable-contraction kernel
    //! (doc/CollocationRewritePlan.md step 5) CANNOT DRIFT APART -- they are the same geometry by
    //! construction, not by two copies of the same twenty lines.  Pure move: every expression and its
    //! evaluation order is unchanged, so \c ForShellPairBox stays bit-identical.
    struct BoxGeom
    {
        bool     live=false;      //!< false = screen (1) killed the whole box; nothing to walk
        double   pMin=0.0;        //!< aMinI+aMinJ, the slowest-decaying product exponent
        double   pfExp=0.0;       //!< the prefactor exponent (screen 1)
        double   reach=0.0;       //!< shrunk exp-tail radius, +1 a.u. polynomial margin
        double   lnCut=0.0;       //!< screen (3) bound
        double   Rq=0.0;          //!< |r-P|^2 chord bound, widened so the bracket is a strict superset
        double   Csz=0.0;         //!< |sz|^2, the chord quadratic's leading coefficient
        rvec3_t  P{0,0,0};        //!< product centre
        rvec3_t  sx{0,0,0}, sy{0,0,0}, sz{0,0,0};   //!< the three grid step vectors
        rvec3_t  r00{0,0,0};      //!< Cartesian position of the box's first corner point
        int      hw[3]={0,0,0};   //!< half-widths in grid indices
        long     c0[3]={0,0,0};   //!< centre grid index
        double   fc[3]={0,0,0};   //!< the product centre in GRID units (fP*N) -- the non-ortho kernel's origin
    };

    //! Build \c BoxGeom for (shell i0, shell j0 imaged to \a Roff) on the \a N-division grid of cell \a A.
    BoxGeom MakeBoxGeom(size_t i0, size_t j0, const rvec3_t& Roff,
                        const UnitCell& A, const ivec3_t& N, double epsEff) const
    {
        BoxGeom g;
        const rvec_t ei=radials[i0]->GetExponents();
        const rvec_t ej=radials[j0]->GetExponents();
        const rvec3_t Ri=radials[i0]->GetCenter(), Rj=radials[j0]->GetCenter()+Roff;
        double aMinI=ei[0]; for (double e:ei) aMinI=std::min(aMinI,e);
        double aMinJ=ej[0]; for (double e:ej) aMinJ=std::min(aMinJ,e);
        g.pMin=aMinI+aMinJ;
        g.P=(aMinI*Ri+aMinJ*Rj)/g.pMin;
        const double lnE=-std::log(epsEff);
        const rvec3_t dij=Ri-Rj;
        g.pfExp=(aMinI*aMinJ/g.pMin)*(dij.x*dij.x+dij.y*dij.y+dij.z*dij.z);
        if (g.pfExp>=lnE) return g;                          // live stays false
        g.reach=std::sqrt((lnE-g.pfExp)/g.pMin)+1.0;
        const rvec3_t fP=A.ToFractional(g.P);
        const rvec3_t fex=A.ToFractional(rvec3_t(1,0,0)), fey=A.ToFractional(rvec3_t(0,1,0)), fez=A.ToFractional(rvec3_t(0,0,1));
        const double rbx=std::sqrt(fex.x*fex.x+fey.x*fey.x+fez.x*fez.x);
        const double rby=std::sqrt(fex.y*fex.y+fey.y*fey.y+fez.y*fez.y);
        const double rbz=std::sqrt(fex.z*fex.z+fey.z*fey.z+fez.z*fez.z);
        g.hw[0]=int(std::ceil(g.reach*rbx*N.x))+1;
        g.hw[1]=int(std::ceil(g.reach*rby*N.y))+1;
        g.hw[2]=int(std::ceil(g.reach*rbz*N.z))+1;
        g.c0[0]=std::lround(fP.x*N.x); g.c0[1]=std::lround(fP.y*N.y); g.c0[2]=std::lround(fP.z*N.z);
        g.fc[0]=fP.x*N.x;              g.fc[1]=fP.y*N.y;              g.fc[2]=fP.z*N.z;
        static const bool kSphere=[]{const char* s=std::getenv("GPW_SPHERE_SCREEN"); return !s || std::atoi(s)!=0;}();
        g.lnCut = kSphere ? std::min(lnE+12.0, g.pMin*g.reach*g.reach+g.pfExp) : lnE+12.0;
        g.sx=A.ToCartesian(rvec3_t(1.0/N.x,0,0));
        g.sy=A.ToCartesian(rvec3_t(0,1.0/N.y,0));
        g.sz=A.ToCartesian(rvec3_t(0,0,1.0/N.z));
        g.r00=A.ToCartesian(rvec3_t(double(g.c0[0]-g.hw[0])/N.x, double(g.c0[1]-g.hw[1])/N.y,
                                    double(g.c0[2]-g.hw[2])/N.z));
        g.Rq=((g.lnCut-g.pfExp)/g.pMin)*(1.0+1e-12)+1e-12;
        g.Csz=g.sz.x*g.sz.x+g.sz.y*g.sz.y+g.sz.z*g.sz.z;
        g.live=true;
        return g;
    }

    // =========================== THE SEPARABLE-CONTRACTION KERNEL ===========================
    // doc/CollocationRewritePlan.md step 5.  Replaces the per-point walk for one (shell pair, offset) with:
    //   1. COLLAPSE the whole shell pair into ONE Gaussian about the product centre times ONE polynomial
    //      (PairPoly) -- the density-matrix weights summed in, so there is no per-component-pair loop left;
    //   2. re-expand that polynomial in GRID-INDEX powers (GridPoly), because on a NON-ORTHO metric -- which
    //      is every production cell (doc there) -- no Cartesian coordinate is a function of one grid index;
    //   3. contract onto the cube with "Mathieu's trick": the quadratic form factors EXACTLY into three 2-D
    //      exponential tables, e^{-p|u|^2} = T12 T23 T31, so the innermost loop is lp+1 fused multiply-adds
    //      and three table lookups -- no exp, no monomial, no component loop.
    // ON BY DEFAULT since 2026-08-27 (GPW_CONTRACT_CUBE=0 opts out); ForShellPairBox remains the reference
    // and src/BasisSet/Molecule/tests/M_PG_BoxWalk.C is the oracle.

    //! \f$(u+p_a)^{n_a}(u+p_b)^{n_b}=\sum_s\alpha_s u^s\f$ -- the ONE-DIMENSIONAL re-expansion about the
    //! product centre.  Per axis, which is exactly what makes the cube a tensor contraction.
    static void Binom1D(int na, int nb, double pa, double pb, double* alpha)
    {
        for (int t=0;t<=na+nb;t++) alpha[t]=0.0;
        for (int k=0;k<=na;k++)
        {
            const double ca=Choose(na,k)*IntPow(pa,na-k);
            for (int l=0;l<=nb;l++) alpha[k+l]+=ca*Choose(nb,l)*IntPow(pb,nb-l);
        }
    }
    static double Choose(int n, int k) { double c=1.0; for (int t=0;t<k;t++) c*=double(n-t)/double(t+1); return c; }
    static double IntPow(double x, int n) { double r=1.0; for (int t=0;t<n;t++) r*=x; return r; }

    static constexpr int kMaxPoly=9;      //!< polynomial extent: lp = La+Lb, so L<=4 per shell
    static constexpr size_t kMaxShell=45; //!< widest Cartesian shell: (L+1)(L+2)/2 = 45 at L=8

    //! The collapsed shell pair: one Gaussian \f$e^{-p|r-P|^2}\f$ times one polynomial in \f$u=r-P\f$,
    //! with the caller's per-component weights already summed in.
    struct PairPoly
    {
        bool    live=false;
        double  p=0.0, Eij=0.0;
        rvec3_t P{0,0,0};
        rvec3_t PA{0,0,0}, PB{0,0,0};   //!< P-A and P-(B+Roff): the binomial shifts, kept so the GATHER
                                        //!< can rebuild the same per-axis coefficients without re-deriving
                                        //!< the image offset (which is not otherwise recoverable from P).
        int     lp=0;
        double  c[kMaxPoly][kMaxPoly][kMaxPoly];
    };

    //! Collapse (shell i0 x shell j0 at \a Roff) weighted by \a w (row-major nI x nJ) into a PairPoly.
    //! \note UNCONTRACTED shells only -- a contracted shell is a SUM over primitive pairs, each with its own
    //! \f$p\f$, \f$P\f$ and box, so it needs one cube per primitive pair (plan open question 2).  Returns
    //! \c live=false otherwise, and the caller falls back to the walk.
    //! \a w may be null: the GATHER direction needs only the Gaussian part (p, P, Eij, lp), since the
    //! per-component coefficients enter there on the way OUT, not on the way in.
    PairPoly MakePairPoly(size_t i0, size_t nI, size_t j0, size_t nJ, const rvec3_t& Roff,
                          const double* w) const
    {
        PairPoly q;
        const rvec_t ei=radials[i0]->GetExponents(), gi=radials[i0]->GetCoeffs();
        const rvec_t ej=radials[j0]->GetExponents(), gj=radials[j0]->GetCoeffs();
        if (ei.size()!=1 || ej.size()!=1) return q;           // contracted: not this kernel
        int lp=0;
        for (size_t a=0;a<nI;a++) lp=std::max(lp, pols[i0+a].GetTotalL());
        int lpj=0;
        for (size_t b=0;b<nJ;b++) lpj=std::max(lpj, pols[j0+b].GetTotalL());
        q.lp=lp+lpj;
        if (q.lp>=kMaxPoly) return q;                         // beyond the tensor extent: fall back
        for (int a=0;a<=q.lp;a++)
            for (int b=0;a+b<=q.lp;b++)
                for (int cc=0;a+b+cc<=q.lp;cc++) q.c[a][b][cc]=0.0;
        // THE GAUSSIAN PRODUCT THEOREM comes from the SHARED GaussProduct, not from a second copy of the
        // same algebra here (user, 2026-08-27).  The contraction coefficients are folded in on top: Omega
        // is a PRIMITIVE-pair object and its Eij carries no g_i g_j.
        const GaussProduct gp(ei[0], radials[i0]->GetCenter(), ej[0], radials[j0]->GetCenter()+Roff);
        q.p=gp.p;
        q.P=gp.P;
        q.Eij=gi[0]*gj[0]*gp.Eij;
        q.PA=gp.PA; q.PB=gp.PB;
        const rvec3_t PA=q.PA, PB=q.PB;
        if (!w) { q.live=true; return q; }                    // Gaussian part only (the gather)
        double ax[kMaxPoly], ay[kMaxPoly], az[kMaxPoly];
        for (size_t a=0;a<nI;a++)
            for (size_t b=0;b<nJ;b++)
            {
                const double wt=w[a*nJ+b]*ns[i0+a]*ns[j0+b];
                if (wt==0.0) continue;
                const Polarization& pi=pols[i0+a];
                const Polarization& pj=pols[j0+b];
                Binom1D(pi.n,pj.n,PA.x,PB.x,ax);
                Binom1D(pi.l,pj.l,PA.y,PB.y,ay);
                Binom1D(pi.m,pj.m,PA.z,PB.z,az);
                for (int sx=0;sx<=pi.n+pj.n;sx++)
                    for (int sy=0;sy<=pi.l+pj.l;sy++)
                    {
                        const double axy=wt*ax[sx]*ay[sy];
                        if (axy==0.0) continue;
                        for (int sz=0;sz<=pi.m+pj.m;sz++) q.c[sx][sy][sz]+=axy*az[sz];
                    }
            }
        q.live=true;
        return q;
    }

    //! A polynomial in the three GRID-INDEX offsets, degree-bounded throughout (a full-extent zero/copy
    //! here costs more than the contraction it feeds -- measured while prototyping).
    struct GridPoly
    {
        double q[kMaxPoly][kMaxPoly][kMaxPoly];
        void Zero(int deg)
        {
            for (int a=0;a<=deg;a++)
                for (int b=0;a+b<=deg;b++)
                    for (int c=0;a+b+c<=deg;c++) q[a][b][c]=0.0;
        }
        void MulLinear(const double L[3], int deg)
        {
            double n[kMaxPoly][kMaxPoly][kMaxPoly];
            for (int a=0;a<=deg+1;a++)
                for (int b=0;a+b<=deg+1;b++)
                    for (int c=0;a+b+c<=deg+1;c++) n[a][b][c]=0.0;
            for (int a=0;a<=deg;a++)
                for (int b=0;a+b<=deg;b++)
                    for (int c=0;a+b+c<=deg;c++)
                    {
                        const double v=q[a][b][c];
                        if (v==0.0) continue;
                        n[a+1][b][c]+=v*L[0]; n[a][b+1][c]+=v*L[1]; n[a][b][c+1]+=v*L[2];
                    }
            for (int a=0;a<=deg+1;a++)
                for (int b=0;a+b<=deg+1;b++)
                    for (int c=0;a+b+c<=deg+1;c++) q[a][b][c]=n[a][b][c];
        }
    };

    //! Re-expand the PairPoly's Cartesian polynomial in powers of the grid-index offsets.
    static GridPoly ToGridPoly(const PairPoly& q, const BoxGeom& g)
    {
        const double Lx[3]={g.sx.x,g.sy.x,g.sz.x};
        const double Ly[3]={g.sx.y,g.sy.y,g.sz.y};
        const double Lz[3]={g.sx.z,g.sy.z,g.sz.z};
        GridPoly out; out.Zero(q.lp);
        GridPoly t;
        for (int sx=0;sx<=q.lp;sx++)
            for (int sy=0;sy+sx<=q.lp;sy++)
                for (int sz=0;sz+sy+sx<=q.lp;sz++)
                {
                    const double cf=q.c[sx][sy][sz];
                    if (cf==0.0) continue;
                    t.Zero(sx+sy+sz); t.q[0][0][0]=cf;
                    int deg=0;
                    for (int r=0;r<sx;r++) { t.MulLinear(Lx,deg); deg++; }
                    for (int r=0;r<sy;r++) { t.MulLinear(Ly,deg); deg++; }
                    for (int r=0;r<sz;r++) { t.MulLinear(Lz,deg); deg++; }
                    for (int a=0;a<=deg;a++)
                        for (int b=0;a+b<=deg;b++)
                            for (int c=0;a+b+c<=deg;c++) out.q[a][b][c]+=t.q[a][b][c];
                }
        return out;
    }

    //! \brief Scatter the collapsed shell pair onto \a dst over its cube.  The NON-ORTHO kernel: three 2-D
    //! exponential tables (Mathieu) and a grid-index polynomial, innermost loop \c lp+1 FMAs.
    //! \brief One row of a 2-D Mathieu table by RECURRENCE: \c out[t] = \f$e^{base+t\,step}\f$.
    //!
    //! CP2K builds these tables this way and we did not (doc/CollocationRewritePlan.md §3c-bis): measured at
    //! the unit level, the three \f$O(n^2)\f$ exponential planes are **84% of the contraction kernel** —
    //! \f$3n^2\f$ scalar `std::exp` calls at ~16 ns each — against 16% for the \f$O(n^3)\f$ point loop.
    //! The exponent is LINEAR in the inner index of each table (only the cross term \f$2e_ae_bh_{ab}\f$
    //! moves), so a row is one seed and \f$n\f$ multiplies: **2 exps per row instead of \f$n\f$.**
    //!
    //! ⚠ SEEDED AT THE LARGEST ENTRY AND WALKED DOWNWARD, ALWAYS.  This is the underflow rule from the
    //! 2026-08-26 analysis, and it is the whole reason the direction is branched on rather than fixed: a
    //! recurrence seeded in the TAIL can start at (or below) the underflow floor and stay identically zero
    //! through entries where the value is significant.  Seeded at the peak, each step only shrinks — an
    //! entry that underflows to 0 is one that `std::exp` would have returned 0 for anyway.
    //! ⚠ And the decrement is \c exp(-step), never \c 1/exp(step): for a steep row \c exp(step) overflows
    //! to infinity and its reciprocal is 0, which would zero the entire row after its seed.
    static void ExpRow(double base, double step, int n, double* out)
    {
        if (step<=0.0)                                   // decreasing in t: t=0 is the largest
        {
            double v=std::exp(base); const double m=std::exp(step);
            for (int t=0;t<n;t++) { out[t]=v; v*=m; }
        }
        else                                             // increasing in t: t=n-1 is the largest
        {
            double v=std::exp(base+double(n-1)*step); const double m=std::exp(-step);
            for (int t=n-1;t>=0;t--) { out[t]=v; v*=m; }
        }
    }
    //! Build the Mathieu tables by \c ExpRow instead of \f$3n^2\f$ \c std::exp calls.  **DEFAULT ON**
    //! (\c GPW_EXP_RECURRENCE=0 opts out); NOT a CP2K deviation — it is what CP2K does.
    //!
    //! ⚠ IT IS NOT BIT-IDENTICAL — a product of \f$n\f$ rounded factors is not the rounded product — so it
    //! was built default-OFF and defaulted ON only on measurement.  What the measurement said:
    //!   - against a NAIVE EXACT reference the contraction goes 1e-15 → **7e-15 relative, FLAT in box size**
    //!     (it is \f$n\varepsilon\f$, and \f$n\f$ is bounded by the box) — still **4× better than the
    //!     WALK's 3e-14**, so the contraction remains the most accurate of the three routes;
    //!   - `ctest -j8` is **793/793 on BOTH settings**;
    //!   - and all three benchmark anchors are unchanged **to all 10 printed s.f.** — Si Γ −7.115067844,
    //!     NaF SR2 Γ −24.4303364755, MnO AFM-II −61.40297551.
    //! ⇒ It is anchor-moving in principle and moved nothing anybody pins.  Keep the knob: it is the A/B for
    //! any future accuracy question, and the first thing to try if a system ever disagrees.
    static bool UseExpRecurrence()
    { static const bool b=[]{const char* s=std::getenv("GPW_EXP_RECURRENCE"); return !s || std::atoi(s)!=0;}(); return b; }

    //! \brief Build the three 2-D Mathieu tables and the \f$e_1\f$ power table for one cube.
    //! ONE definition, shared by \c ContractCubeN and \c GatherCubeN — they held two verbatim copies, and a
    //! second copy of a rule is the defect class that has already cost this tree twice (the integrate-back
    //! \c Re[] screen, the memo logic).  The two directions MUST see identical tables or they stop being
    //! adjoints, which is now true by construction rather than by inspection.
    void BuildCubeTables(const BoxGeom& g, double p, int lp, const int n[3], const double h[3][3],
                         rvec_t& T12, rvec_t& T23, rvec_t& T31, rvec_t& E1) const
    {
        T12.resize(size_t(n[0])*n[1]); T23.resize(size_t(n[1])*n[2]); T31.resize(size_t(n[2])*n[0]);
        E1.resize(size_t(lp+1)*n[0]);
        auto eOf=[&](int a, int t){ return double(g.c0[a]-g.hw[a]+t)-g.fc[a]; };
        if (UseExpRecurrence())
        {
            // The inner index steps e_b by exactly 1, so the exponent steps by the constant -2 p e_a h_ab.
            for (int i=0;i<n[0];i++)
            { const double e1=eOf(0,i);
              ExpRow(-p*(e1*e1*h[0][0]+2*e1*eOf(1,0)*h[0][1]), -2.0*p*e1*h[0][1], n[1], &T12[size_t(i)*n[1]]); }
            for (int j=0;j<n[1];j++)
            { const double e2=eOf(1,j);
              ExpRow(-p*(e2*e2*h[1][1]+2*e2*eOf(2,0)*h[1][2]), -2.0*p*e2*h[1][2], n[2], &T23[size_t(j)*n[2]]); }
            for (int k=0;k<n[2];k++)
            { const double e3=eOf(2,k);
              ExpRow(-p*(e3*e3*h[2][2]+2*e3*eOf(0,0)*h[2][0]), -2.0*p*e3*h[2][0], n[0], &T31[size_t(k)*n[0]]); }
        }
        else
        {
            for (int i=0;i<n[0];i++) for (int j=0;j<n[1];j++)
            { const double e1=eOf(0,i), e2=eOf(1,j); T12[size_t(i)*n[1]+j]=std::exp(-p*(e1*e1*h[0][0]+2*e1*e2*h[0][1])); }
            for (int j=0;j<n[1];j++) for (int k=0;k<n[2];k++)
            { const double e2=eOf(1,j), e3=eOf(2,k); T23[size_t(j)*n[2]+k]=std::exp(-p*(e2*e2*h[1][1]+2*e2*e3*h[1][2])); }
            for (int k=0;k<n[2];k++) for (int i=0;i<n[0];i++)
            { const double e3=eOf(2,k), e1=eOf(0,i); T31[size_t(k)*n[0]+i]=std::exp(-p*(e3*e3*h[2][2]+2*e3*e1*h[2][0])); }
        }
        for (int i=0;i<n[0];i++) { double e=1.0; const double e1=eOf(0,i);
                                   for (int a=0;a<=lp;a++) { E1[size_t(a)*n[0]+i]=e; e*=e1; } }
    }
    //! \brief The \f$l_p\f$-SPECIALIZED body.  \a LP is \c q.lp, and every loop bound that was a runtime
    //! `q.lp` is now a compile-time constant: the innermost accumulation unrolls, \c ci[] can live in
    //! registers across the point loop, and the fold arrays shrink from \c kMaxPoly to what is used.
    //! ⚠ BIT-IDENTICAL BY CONSTRUCTION — a pure substitution of one constant for one variable.  The terms
    //! and their summation ORDER are untouched (C++ may not reassociate floating point), which is what lets
    //! this land outside an anchor re-bank.  Dropping the \c E1 table for Horner would NOT be, and is a
    //! separate increment (doc/CollocationRewritePlan.md §3c-bis).
    template <int LP>
    void ContractCubeN(const BoxGeom& g, const PairPoly& q, const ivec3_t& N, rvec_t& dst) const
    {
        const int n[3]={2*g.hw[0]+1, 2*g.hw[1]+1, 2*g.hw[2]+1};
        const GridPoly Q=ToGridPoly(q,g);
        // the metric h_ab = s_a . s_b, and e_a = index_a - fc[a]
        const rvec3_t sv[3]={g.sx,g.sy,g.sz};
        double h[3][3];
        for (int a=0;a<3;a++)
            for (int b=0;b<3;b++) h[a][b]=sv[a].x*sv[b].x+sv[a].y*sv[b].y+sv[a].z*sv[b].z;
        auto eOf=[&](int a, int t){ return double(g.c0[a]-g.hw[a]+t)-g.fc[a]; };
        static thread_local rvec_t T12,T23,T31,E1;
        BuildCubeTables(g,q.p,LP,n,h,T12,T23,T31,E1);            // ONE definition -- see BuildCubeTables
        // ⚠ THE O(n^2) COST IS THE PER-LINE WORK, NOT THE TABLES -- measured, and it refuted a fit
        // (M_PG_BoxWalk.WhereTheContractionSpendsItsTime, 2026-08-28).  A raster-scaling fit attributed 56%
        // of the kernel to an O(N^2) term and the obvious suspect was the three exponential tables; calling
        // BuildCubeTables directly says they are **7%**.  What is really O(n^2) is THIS: the (j,k) line
        // setup, run n1*n2 ~ 1849 times per task -- the e2 fold, a sqrt, and the modulo wraps.  Two cheap
        // consequences, both bit-identical:
        //   (1) WRAP TABLES.  The wrapped raster index of an axis depends only on that axis's index, so it
        //       is O(n) to tabulate and O(n^2) to recompute.  The per-line integer modulos (a hardware
        //       divide each -- the same instruction the 2026-08-26 box-walk work took off the inner path)
        //       become loads from an n-element table built once per task.
        //   (2) CHORD BEFORE FOLD.  ~21% of (j,k) lines miss the sphere entirely (pi/4 of a square), and
        //       the e2 fold used to run BEFORE the test that discovers it.  Test first, fold only if the
        //       line has points.
        static thread_local std::vector<size_t> wrapX, wrapY, wrapZ;
        auto wrapAxis=[&](int a, long N_a, std::vector<size_t>& w)
        {
            w.resize(size_t(n[a]));
            for (int t=0;t<n[a];t++) w[size_t(t)]=size_t((((g.c0[a]-g.hw[a]+t)%N_a)+N_a)%N_a);
        };
        wrapAxis(0,N.x,wrapX); wrapAxis(1,N.y,wrapY); wrapAxis(2,N.z,wrapZ);
        double cij[LP+1][LP+1], ci[LP+1];
                // ⚠ HOISTED, and the reason is worth stating (user, 2026-08-28: *"the quadratic looks
                // like something that does not depend on rho(r), and therefore does not change with
                // iterations"*).  Almost right: h_ab and the e's ARE pure geometry, and the ONE
                // D-dependent term is Rq, whose eps comes from the D-aware tolerance -- so the sphere's
                // RADIUS breathes with the density and nothing else does.  Caching the geometry half is
                // O(n^2) PER TASK, i.e. the value-cache trade step 7 deleted 4 GB of.  What IS free is
                // this: several of these terms do not vary with j at all, and were being recomputed on
                // every one of ~1849 lines -- including a floating-point DIVIDE for a task-invariant
                // constant.  All three hoists are bit-identical (same expressions, same association,
                // evaluated once instead of many times).
        const double inv=0.5/h[0][0];                      // task-invariant: a DIVIDE, not a per-line one
        for (int k=0;k<n[2];k++)
        {
            const double e3=eOf(2,k);
            const double b3=h[2][0]*e3;                    // per-PLANE halves of B2 and C2
            const double c3=h[2][2]*e3*e3;
            for (int a=0;a<=LP;a++)                        // fold e3 once per plane
                for (int b=0;a+b<=LP;b++)
                {
                    double v=0.0, e=1.0;
                    for (int c=0;a+b+c<=LP;c++) { v+=Q.q[a][b][c]*e; e*=e3; }
                    cij[a][b]=v;
                }
            const size_t mk=wrapZ[size_t(k)];
            for (int j=0;j<n[1];j++)
            {
                const double e2=eOf(1,j);
                // the chord in e1: h00 e1^2 + 2 e1 (h01 e2 + h20 e3) + (h11 e2^2 + h22 e3^2 + 2 h12 e2 e3) <= Rq
                const double B2=2.0*(h[0][1]*e2+b3);
                const double C2=h[1][1]*e2*e2+c3+2.0*h[1][2]*e2*e3-g.Rq;
                const double D2=B2*B2-4.0*h[0][0]*C2;
                if (D2<0.0) continue;                        // the line misses the sphere
                const double sq=std::sqrt(D2);
                // THE EXACT CHORD.  A point t is kept iff the quadratic is <=0 there, i.e. iff
                // ceil(root1) <= t <= floor(root2) -- so those ARE the bounds.  A first cut used
                // floor(root1)-1 .. ceil(root2)+1, which keeps up to three extra points at each end; for
                // one pair those all carry the SAME SIGN, so they do not cancel in an integral.  Measured:
                // the integrated per-pair difference was ~30x the pointwise one and accumulated over ~400
                // pairs into 6.2e-6 Ha on NaF.  The +-kChordFuzz guard (1e-9 of an index) keeps a point
                // that sits on the boundary rather than dropping it to a rounding difference.
                constexpr double kChordFuzz=1e-9;
                long ia=long(std::ceil ((-B2-sq)*inv+g.fc[0]-kChordFuzz))-(g.c0[0]-g.hw[0]);
                long ib=long(std::floor((-B2+sq)*inv+g.fc[0]+kChordFuzz))-(g.c0[0]-g.hw[0]);
                if (ia<0) ia=0;
                if (ib>n[0]-1) ib=n[0]-1;
                if (ia>ib) continue;
                for (int a=0;a<=LP;a++)                    // fold e2 -- ONLY for a line with points
                {
                    double v=0.0, e=1.0;
                    for (int b=0;a+b<=LP;b++) { v+=cij[a][b]*e; e*=e2; }
                    ci[a]=v;
                }
                const double t23=T23[size_t(j)*n[2]+k];
                const size_t my=wrapY[size_t(j)];
                size_t mi=wrapX[size_t(ia)];
                for (long ii=ia; ii<=ib; ii++, mi=(mi+1==size_t(N.x)?0:mi+1))
                {
                    double v=0.0;
                    for (int a=0;a<=LP;a++) v+=ci[a]*E1[size_t(a)*n[0]+ii];   // <- lp+1 FMAs
                    dst[(mi*size_t(N.y)+my)*size_t(N.z)+mk]
                        += q.Eij*v*T12[size_t(ii)*n[1]+j]*t23*T31[size_t(k)*n[0]+ii];
                }
            }
        }
    }
    //! \brief The GATHER: the exact TRANSPOSE of ContractCube.  Accumulates the grid-index MOMENTS
    //! \f$W_{abc}=\sum_g V(g)\,e_1^ae_2^be_3^c\,e^{-p|u|^2}\f$ over the same cube, same tables, same chord.
    //! Because it is the transpose of the scatter -- not merely "the same physics backwards" -- the
    //! collocate/integrate pair stays an EXACT ADJOINT, which is what makes the GPW energy variational.
    //! The \f$l_p\f$-specialized gather body; see \c ContractCubeN for why it is bit-identical.
    template <int LP>
    void GatherCubeN(const BoxGeom& g, const PairPoly& q, const ivec3_t& N, const rvec_t& V,
                     GridPoly& W) const
    {
        const int n[3]={2*g.hw[0]+1, 2*g.hw[1]+1, 2*g.hw[2]+1};
        const rvec3_t sv[3]={g.sx,g.sy,g.sz};
        double h[3][3];
        for (int a=0;a<3;a++)
            for (int b=0;b<3;b++) h[a][b]=sv[a].x*sv[b].x+sv[a].y*sv[b].y+sv[a].z*sv[b].z;
        auto eOf=[&](int a, int t){ return double(g.c0[a]-g.hw[a]+t)-g.fc[a]; };
        static thread_local rvec_t T12,T23,T31,E1;
        BuildCubeTables(g,q.p,LP,n,h,T12,T23,T31,E1);            // ONE definition -- see BuildCubeTables
        W.Zero(LP);
        const double inv=0.5/h[0][0];                      // hoisted -- see ContractCubeN
        // Wrap tables, exactly as in ContractCubeN (see the note there): the per-line integer modulos are
        // O(n^2) hardware divides for an O(n) fact.  The gather already tests the chord before its fold.
        static thread_local std::vector<size_t> wrapX, wrapY, wrapZ;
        auto wrapAxis=[&](int a, long N_a, std::vector<size_t>& w)
        {
            w.resize(size_t(n[a]));
            for (int t=0;t<n[a];t++) w[size_t(t)]=size_t((((g.c0[a]-g.hw[a]+t)%N_a)+N_a)%N_a);
        };
        wrapAxis(0,N.x,wrapX); wrapAxis(1,N.y,wrapY); wrapAxis(2,N.z,wrapZ);
        double cij[LP+1][LP+1], ci[LP+1];
        for (int k=0;k<n[2];k++)
        {
            const double e3=eOf(2,k);
            for (int a=0;a<=LP;a++) for (int b=0;a+b<=LP;b++) cij[a][b]=0.0;
            const size_t mk=wrapZ[size_t(k)];
            const double b3=h[2][0]*e3, c3=h[2][2]*e3*e3;   // per-PLANE halves of B2 and C2
            for (int j=0;j<n[1];j++)
            {
                const double e2=eOf(1,j);
                const double B2=2.0*(h[0][1]*e2+b3);
                const double C2=h[1][1]*e2*e2+c3+2.0*h[1][2]*e2*e3-g.Rq;
                const double D2=B2*B2-4.0*h[0][0]*C2;
                if (D2<0.0) continue;
                const double sq=std::sqrt(D2);
                // THE EXACT CHORD.  A point t is kept iff the quadratic is <=0 there, i.e. iff
                // ceil(root1) <= t <= floor(root2) -- so those ARE the bounds.  A first cut used
                // floor(root1)-1 .. ceil(root2)+1, which keeps up to three extra points at each end; for
                // one pair those all carry the SAME SIGN, so they do not cancel in an integral.  Measured:
                // the integrated per-pair difference was ~30x the pointwise one and accumulated over ~400
                // pairs into 6.2e-6 Ha on NaF.  The +-kChordFuzz guard (1e-9 of an index) keeps a point
                // that sits on the boundary rather than dropping it to a rounding difference.
                constexpr double kChordFuzz=1e-9;
                long ia=long(std::ceil ((-B2-sq)*inv+g.fc[0]-kChordFuzz))-(g.c0[0]-g.hw[0]);
                long ib=long(std::floor((-B2+sq)*inv+g.fc[0]+kChordFuzz))-(g.c0[0]-g.hw[0]);
                if (ia<0) ia=0;
                if (ib>n[0]-1) ib=n[0]-1;
                if (ia>ib) continue;
                for (int a=0;a<=LP;a++) ci[a]=0.0;
                const double t23=T23[size_t(j)*n[2]+k];
                const size_t my=wrapY[size_t(j)];
                size_t mi=wrapX[size_t(ia)];
                for (long ii=ia; ii<=ib; ii++, mi=(mi+1==size_t(N.x)?0:mi+1))
                {
                    const double gv=V[(mi*size_t(N.y)+my)*size_t(N.z)+mk]
                                    *T12[size_t(ii)*n[1]+j]*t23*T31[size_t(k)*n[0]+ii];
                    for (int a=0;a<=LP;a++) ci[a]+=gv*E1[size_t(a)*n[0]+ii];   // <- lp+1 FMAs
                }
                for (int a=0;a<=LP;a++)                    // unfold e2, once per line
                {
                    double e=1.0;
                    for (int b=0;a+b<=LP;b++) { cij[a][b]+=ci[a]*e; e*=e2; }
                }
            }
            for (int a=0;a<=LP;a++)                        // unfold e3, once per plane
                for (int b=0;a+b<=LP;b++)
                {
                    double e=1.0;
                    for (int c=0;a+b+c<=LP;c++) { W.q[a][b][c]+=cij[a][b]*e; e*=e3; }
                }
        }
    }

    //! \brief THE \f$l_p\f$ DISPATCH (doc/CollocationRewritePlan.md §3c-bis, stage 1).  One arm per
    //! degree, so the grid loops run with \f$l_p\f$ known at compile time — CP2K dispatches its collocation
    //! the same way (`grid_cpu_collocate.c`, `always_inline` low-\f$l_p\f$ variants).
    //!
    //! ⚠ WHAT THIS CAN AND CANNOT BUY, measured before it was built
    //! (`M_PG_BoxWalk.WhereTheContractionSpendsItsTime`): with the box held fixed,
    //! \f$\text{contract}(l_p)=116.7+9.7(l_p{+}1)\f$ µs, so the \f$l_p\f$-DEPENDENT work this
    //! specializes is **20% of the kernel at \f$l_p{=}2\f$ and 29% at \f$l_p{=}4\f$**.  The other
    //! ~63% is the three 2-D exponential tables — \f$3n^2\f$ scalar `std::exp` calls — and no amount of
    //! unrolling touches them.  Do not quote this dispatch as "the collocation got faster"; quote the two
    //! ledger buckets.
    //!
    //! The switch is TOTAL: \c MakePairPoly declines \f$l_p\ge\f$ \c kMaxPoly and returns \c live=false,
    //! so a caller that reached here has \f$0\le l_p<9\f$.
    void ContractCube(const BoxGeom& g, const PairPoly& q, const ivec3_t& N, rvec_t& dst) const
    {
        switch (q.lp)
        {
            case 0: ContractCubeN<0>(g,q,N,dst); return;
            case 1: ContractCubeN<1>(g,q,N,dst); return;
            case 2: ContractCubeN<2>(g,q,N,dst); return;
            case 3: ContractCubeN<3>(g,q,N,dst); return;
            case 4: ContractCubeN<4>(g,q,N,dst); return;
            case 5: ContractCubeN<5>(g,q,N,dst); return;
            case 6: ContractCubeN<6>(g,q,N,dst); return;
            case 7: ContractCubeN<7>(g,q,N,dst); return;
            case 8: ContractCubeN<8>(g,q,N,dst); return;
        }
        assert(false && "lp outside [0,kMaxPoly) -- MakePairPoly should have declined this pair");
    }
    //! The gather's dispatch; the two MUST stay in step or they stop being adjoints.
    void GatherCube(const BoxGeom& g, const PairPoly& q, const ivec3_t& N, const rvec_t& V,
                    GridPoly& W) const
    {
        switch (q.lp)
        {
            case 0: GatherCubeN<0>(g,q,N,V,W); return;
            case 1: GatherCubeN<1>(g,q,N,V,W); return;
            case 2: GatherCubeN<2>(g,q,N,V,W); return;
            case 3: GatherCubeN<3>(g,q,N,V,W); return;
            case 4: GatherCubeN<4>(g,q,N,V,W); return;
            case 5: GatherCubeN<5>(g,q,N,V,W); return;
            case 6: GatherCubeN<6>(g,q,N,V,W); return;
            case 7: GatherCubeN<7>(g,q,N,V,W); return;
            case 8: GatherCubeN<8>(g,q,N,V,W); return;
        }
        assert(false && "lp outside [0,kMaxPoly) -- MakePairPoly should have declined this pair");
    }

    //! \brief Turn the grid-index moments \a W into the per-component-pair integrals
    //! \f$b_{ab}=\sum_g V\chi_i\chi_j\f$ (row-major nI x nJ into \a bOut).  This is the transpose of the
    //! collapse: each Cartesian monomial's grid-index polynomial is contracted against \a W, then the same
    //! per-axis binomials that BUILT the coefficient tensor read it back out.
    void MomentsToPairs(const BoxGeom& g, const PairPoly& q, const GridPoly& W,
                        size_t i0, size_t nI, size_t j0, size_t nJ, double* bOut) const
    {
        const double Lx[3]={g.sx.x,g.sy.x,g.sz.x};
        const double Ly[3]={g.sx.y,g.sy.y,g.sz.y};
        const double Lz[3]={g.sx.z,g.sy.z,g.sz.z};
        // M[sx][sy][sz] = Eij * <grid-index polynomial of u^s , W>
        static thread_local double M[kMaxPoly][kMaxPoly][kMaxPoly];
        GridPoly t;
        for (int sx=0;sx<=q.lp;sx++)
            for (int sy=0;sy+sx<=q.lp;sy++)
                for (int sz=0;sz+sy+sx<=q.lp;sz++)
                {
                    t.Zero(sx+sy+sz); t.q[0][0][0]=1.0;
                    int deg=0;
                    for (int r=0;r<sx;r++) { t.MulLinear(Lx,deg); deg++; }
                    for (int r=0;r<sy;r++) { t.MulLinear(Ly,deg); deg++; }
                    for (int r=0;r<sz;r++) { t.MulLinear(Lz,deg); deg++; }
                    double v=0.0;
                    for (int a=0;a<=deg;a++)
                        for (int b=0;a+b<=deg;b++)
                            for (int c=0;a+b+c<=deg;c++) v+=t.q[a][b][c]*W.q[a][b][c];
                    M[sx][sy][sz]=q.Eij*v;
                }
        const rvec3_t PA=q.PA, PB=q.PB;                       // carried on the PairPoly, image offset included
        double ax[kMaxPoly], ay[kMaxPoly], az[kMaxPoly];
        for (size_t a=0;a<nI;a++)
            for (size_t b=0;b<nJ;b++)
            {
                const Polarization& pi=pols[i0+a];
                const Polarization& pj=pols[j0+b];
                Binom1D(pi.n,pj.n,PA.x,PB.x,ax);
                Binom1D(pi.l,pj.l,PA.y,PB.y,ay);
                Binom1D(pi.m,pj.m,PA.z,PB.z,az);
                double v=0.0;
                for (int sx=0;sx<=pi.n+pj.n;sx++)
                    for (int sy=0;sy<=pi.l+pj.l;sy++)
                    {
                        const double axy=ax[sx]*ay[sy];
                        if (axy==0.0) continue;
                        for (int sz=0;sz<=pi.m+pj.m;sz++) v+=axy*az[sz]*M[sx][sy][sz];
                    }
                bOut[a*nJ+b]=ns[i0+a]*ns[j0+b]*v;
            }
    }

    //! The COLLOCATION KERNEL, both directions at once (they must agree or the adjoint breaks).  DEFAULT ON
    //! since 2026-08-27: the separable contraction is 3.6-5.0x faster than the walk AND measurably more
    //! accurate (it carries no per-component screen -- doc/CollocationRewritePlan.md steps 5+6), and step 7
    //! deleted the value cache that used to hide the walk's cost, so the walk is no longer an affordable
    //! default.  \c GPW_CONTRACT_CUBE=0 is the opt-OUT, back onto \c ForShellPairBox -- kept because that
    //! walk is the REFERENCE implementation the unit oracle in \c M_PG_BoxWalk.C checks the kernel against.
    static bool UseContractCube()
    { static const bool b=[]{const char* s=std::getenv("GPW_CONTRACT_CUBE"); return !s || std::atoi(s)!=0;}(); return b; }

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
        // The box geometry now lives in MakeBoxGeom, so the separable-contraction kernel shares it rather
        // than re-deriving it (see BoxGeom).  Identical expressions in identical order -> bit-identical.
        const BoxGeom bg=MakeBoxGeom(i0,j0,Roff,A,N,epsEff);
        if (!bg.live) return;                                // screen (1): the whole box is below eps
        const double pMin=bg.pMin;
        const rvec3_t P=bg.P;
        const double reach=bg.reach;
        const int hwx=bg.hw[0], hwy=bg.hw[1], hwz=bg.hw[2];
        const long cx=bg.c0[0], cy=bg.c0[1], cz=bg.c0[2];
        // Screen (3) is the REACH SPHERE the box is already the BOUNDING BOX OF -- capped by the historical
        // fixed log margin lnE+12 (2026-08-26).  The old bound was lnE+12 alone, and for a DIFFUSE pair that
        // exceeds the box entirely (pMin=0.3 gives r_screen=10.8 au against reach=9.76): the screen never
        // fired, so the CORNERS of the bounding box -- ~48% of its points, every one with an exponential
        // envelope already below epsEff -- paid the full exp/poly evaluation, and the consumer's own
        // |val|<eps test then threw the results away.  Screening on the sphere makes the kept set ISOTROPIC
        // instead of a rectangular hull (which was an artifact of the walk order, not of any tolerance).
        //
        // THE MARGIN IS NOT LOST, it is the one the box already grants: pMin*reach^2+pf expands to
        // lnE + 2 sqrt(pMin (lnE-pf)) + pMin, i.e. reach's +1 a.u. polynomial margin re-expressed in logs
        // (+5.5 at pMin=0.3), and the min() keeps the old +12 wherever it is the tighter of the two (sharp
        // pairs).  So this drops ONLY terms the axis-wise box truncation was already dropping at the same
        // radius.  Measured: Si and NaF total energies unchanged to every printed digit on BOTH the cached
        // and the uncached path (NaF to 12 s.f.), against a cached-vs-uncached spread of 1.25e-6 that this
        // sits far beneath; 771/771.  GPW_SPHERE_SCREEN=0 restores the rectangular hull for A/B.
        const double lnCut=bg.lnCut;                         // screen (3), see BoxGeom
        const rvec3_t sx=bg.sx, sy=bg.sy, sz=bg.sz;          // the constant per-index lattice steps
        const rvec3_t r00=bg.r00;                            // the box's first corner
        const double Rq=bg.Rq, Csz=bg.Csz;                   // the chord bound and its leading coefficient
        // The per-component factor arrays, allocated ONCE per box (kMaxShell bounds a Cartesian shell:
        // (L+1)(L+2)/2 = 45 at L=8, well past any basis this evaluator sees).
        assert(nI<=kMaxShell && nJ<=kMaxShell && "ForShellPairBox: shell wider than kMaxShell components");
        double fI[kMaxShell], fJ[kMaxShell];
        // MONOMIAL POWER TABLES (2026-08-26).  The point body used to evaluate pols[i](d) per component --
        // three intpow/uintpow calls each -- and BOTH Polarization::operator() and uintpow are out-of-line
        // CROSS-DSO calls from here (perf on the uncached Si kernel: Polarization 12.1% + 1.7% @plt,
        // uintpow 6.1% + 1.2% @plt = ~21% of the walk spent producing x^n y^l z^m).  But a shell's
        // components are the monomials of ONE degree over the SAME displacement, so the per-axis powers
        // x^0..x^{Lx} are SHARED across them: build the three tables once per point and index them.
        // PowerTable reproduces uintpow's exact associations, so every fI/fJ is BIT-IDENTICAL.
        static constexpr int kMaxDeg=8;                     // (L+1)(L+2)/2 = 45 = kMaxShell at L=8
        int pnI[kMaxShell], plI[kMaxShell], pmI[kMaxShell], pnJ[kMaxShell], plJ[kMaxShell], pmJ[kMaxShell];
        int LxI=0,LyI=0,LzI=0, LxJ=0,LyJ=0,LzJ=0;
        for (size_t a=0;a<nI;a++)
        {
            const Polarization& p=pols[i0+a];
            pnI[a]=p.n; plI[a]=p.l; pmI[a]=p.m;
            LxI=std::max(LxI,p.n); LyI=std::max(LyI,p.l); LzI=std::max(LzI,p.m);
        }
        for (size_t b=0;b<nJ;b++)
        {
            const Polarization& p=pols[j0+b];
            pnJ[b]=p.n; plJ[b]=p.l; pmJ[b]=p.m;
            LxJ=std::max(LxJ,p.n); LyJ=std::max(LyJ,p.l); LzJ=std::max(LzJ,p.m);
        }
        assert(LxI<=kMaxDeg && LyI<=kMaxDeg && LzI<=kMaxDeg && LxJ<=kMaxDeg && LyJ<=kMaxDeg && LzJ<=kMaxDeg
               && "ForShellPairBox: monomial degree beyond kMaxDeg");
        double pxI[kMaxDeg+1], pyI[kMaxDeg+1], pzI[kMaxDeg+1], pxJ[kMaxDeg+1], pyJ[kMaxDeg+1], pzJ[kMaxDeg+1];
        // INCREMENTAL MODULO WRAP (2026-08-26).  The wrap used to cost SIX integer divisions per point
        // (((g%N)+N)%N on each axis) even though the grid index advances by exactly ONE per loop step.
        // Starting from an in-range residue and bumping it by 1 with a compare reproduces the same
        // residue for ANY number of wraps, so mx/my hoist to their own loops (their residue is constant
        // across the inner ones) and the row offset (mx*N.y+my)*N.z hoists with them -- the innermost
        // body keeps one compare.  Integer arithmetic: identical indices, not merely equivalent ones.
        auto wrap0=[](long g, int n){ return size_t(((g%n)+n)%n); };   // once per axis per loop, not per point
        rvec3_t rx=r00;
        size_t mx=wrap0(cx-hwx,N.x);
        for (int dx=-hwx; dx<=hwx; dx++, rx=rx+sx, mx=(mx+1==size_t(N.x)?0:mx+1))
        {
            rvec3_t ry=rx;
            size_t my=wrap0(cy-hwy,N.y);
            const size_t planeBase=mx*size_t(N.y);
            for (int dy=-hwy; dy<=hwy; dy++, ry=ry+sy, my=(my+1==size_t(N.y)?0:my+1))
            {
                // Solve this line's chord: t^2 Csz + t bq + cq <= 0, t = dz+hwz in [0, 2hwz].
                const rvec3_t g0=ry-P;
                const double bq=2.0*(g0.x*sz.x+g0.y*sz.y+g0.z*sz.z);
                const double cq=(g0.x*g0.x+g0.y*g0.y+g0.z*g0.z)-Rq;
                const double Dq=bq*bq-4.0*Csz*cq;
                if (Dq<0.0) continue;                                  // the line misses the sphere entirely
                const double sq=std::sqrt(Dq), inv=0.5/Csz;
                int t0=int(std::floor((-bq-sq)*inv))-1;                // widened by a point at each end
                int t1=int(std::ceil ((-bq+sq)*inv))+1;
                if (t0<0) t0=0;
                if (t1>2*hwz) t1=2*hwz;
                if (t0>t1) continue;
                rvec3_t r=ry;
                size_t mz=wrap0(cz-hwz,N.z);
                const size_t rowBase=(planeBase+my)*size_t(N.z);
                // Walk the HEAD the cheap way.  r is still ACCUMULATED point by point (never jumped to
                // ry + t0*sz, which would round differently and, being a z-only shortcut, would make the
                // walk anisotropic -- the defect that sank the exp recurrence, doc/OpenWork.md).  So every
                // point the body sees carries exactly the r it carried before.
                for (int t=0;t<t0;t++) { r=r+sz; mz=(mz+1==size_t(N.z)?0:mz+1); }
                for (int dz=-hwz+t0; dz<=-hwz+t1; dz++, r=r+sz, mz=(mz+1==size_t(N.z)?0:mz+1))
                {
                    const rvec3_t di=r-Ri, dj=r-Rj;
                    const double ri2=di.x*di.x+di.y*di.y+di.z*di.z, rj2=dj.x*dj.x+dj.y*dj.y+dj.z*dj.z;
                    if (aMinI*ri2+aMinJ*rj2 > lnCut) continue;   // slowest term < eps*e^-12 -> skip exp/poly evals
                    double radI=0.0; for (size_t p=0;p<ei.size();p++) radI+=gi[p]*std::exp(-ei[p]*ri2);
                    double radJ=0.0; for (size_t q=0;q<ej.size();q++) radJ+=gj[q]*std::exp(-ej[q]*rj2);
                    // Associated exactly as the unblocked (ni*pols[i](di)*radI)*(nj*pols[j](dj)*radJ) was.
                    PowerTable(di.x,LxI,pxI); PowerTable(di.y,LyI,pyI); PowerTable(di.z,LzI,pzI);
                    PowerTable(dj.x,LxJ,pxJ); PowerTable(dj.y,LyJ,pyJ); PowerTable(dj.z,LzJ,pzJ);
                    for (size_t a=0;a<nI;a++) fI[a]=ns[i0+a]*((pxI[pnI[a]]*pyI[plI[a]])*pzI[pmI[a]])*radI;
                    for (size_t b=0;b<nJ;b++) fJ[b]=ns[j0+b]*((pxJ[pnJ[b]]*pyJ[plJ[b]])*pzJ[pmJ[b]])*radJ;
                    assert(mz==size_t((((cz+dz)%N.z)+N.z)%N.z) && "incremental z wrap drifted from the modulo");
                    f(rowBase+mz, fI, fJ);
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
    // --- THE BOX TASK LIST for the analytic collocation / integrate-back -----------------------------------
    // doc/CollocationRewritePlan.md 3c.  A TASK is one (shell pair, cross-cell offset): the unit of work that
    // the box walk and the separable-contraction kernel both consume.  Enumerating and screening them is pure
    // GEOMETRY -- identical across SCF iterations, across k-blocks and across the two directions -- so it is
    // derived ONCE per basis instance (the cell and the centres ARE this basis's own data) and only the
    // ARITHMETIC repeats.
    //
    // THIS REPLACED A 4 GB VALUE CACHE (2026-08-27, plan step 7).  The old design stored the per-point pair
    // VALUES, O(n^3) per (pair, offset) -- 3.9 GB on the MnO magnetic cell -- because the on-the-fly walk was
    // ~100x off CP2K's and re-evaluating every iteration was unaffordable.  Once the separable contraction
    // landed the same 3.9 GB bought only ~1.1x on the run (2.9x on the two box-walk buckets, which are now
    // ~58% of it), so the trade stopped being defensible: measured on the MnO acceptance probe, peak RSS
    // 3915 -> 155 MB against ~1.1x CPU.  The task list is the O(1)-per-task residue of that cache -- ~40 B
    // against ~100 kB per task -- and it is what makes "re-evaluate every iteration" a DESIGN rather than a
    // regression.
    //
    // THE LIST IS NOT WHERE THE TIME WAS, and that is measured, not assumed
    // (M_PG_BoxWalk.OffsetEnumerationIsNotTheCost): enumerating + screening all 4482 tasks of an MnO-shaped
    // fixture costs 2.26 ms against 2294 ms to walk them -- 0.10% of the pass.  So this hoists the geometry
    // for the reasons above and for the longest-first parallel order below, NOT to save the enumeration.
    //
    // D NEVER ENTERS IT (plan 2a).  The list is screened at the bare tolerance kDensityEps; a per-iteration
    // density weight can only RAISE a pair's tolerance (eps_ij = eps/|c_ij| >= eps), hence only ever REMOVE
    // tasks, never add one.  The static list is therefore a strict SUPERSET, and the D-aware screen stays an
    // O(1) per-task predicate outside the walk.  The T3 fold is likewise a per-iteration filter, so the list
    // is fold-independent and survives every fold change.
    struct BoxTask
    {
        ivec3_t n{0,0,0};      //!< the integer cell offset (the caller's Bloch-phase key)
        rvec3_t Roff{0,0,0};   //!< its Cartesian translation
        double  pf=0.0;        //!< PairPrefactorExp: screen (1), a shell-pair property of this offset
    };
    struct ShellPairTasks { std::vector<BoxTask> tasks; };
    //! The tasks of every shell pair, indexed exactly as the (a<=b) shell-pair loops enumerate them.  Built
    //! lazily on first use and never invalidated -- \a A and the centres are this basis instance's own data,
    //! the same assumption the integrate-back memo makes.
    const std::vector<ShellPairTasks>& BoxTasks(const UnitCell& A) const
    {
        if (!itsBoxTasks.empty()) return itsBoxTasks;
        qchem::report::Timed timer("setup: box task list");
        const std::vector<Shell>& shells=Shells();
        const double lnE=-std::log(kDensityEps());
        size_t nTasks=0;
        for (size_t a=0;a<shells.size();a++)
            for (size_t b=a;b<shells.size();b++)
            {
                itsBoxTasks.emplace_back();
                std::vector<BoxTask>& t=itsBoxTasks.back().tasks;
                ForImageOffsets(shells[a].begin,shells[b].begin,A,[&](const ivec3_t& n, const rvec3_t& Roff)
                {
                    const double pf=PairPrefactorExp(shells[a].begin,shells[b].begin,Roff);
                    if (pf<lnE) t.push_back({n,Roff,pf});   // else the whole box is sub-eps for EVERY D
                });
                nTasks+=t.size();
            }
        // PROVENANCE, unconditionally (doc/Benchmark.md's standing rule): a timing row is not reproducible
        // unless it says WHICH collocation kernel produced it.  Printed from the kernel's own side rather
        // than from the run banner, so there is no second copy of the default rule to drift (it is NOT a
        // CP2K deviation -- CP2K collocates exactly this way -- so it does not belong on that table).
        std::cout<<"[collocation] kernel="<<(UseContractCube() ? "separable contraction"
                                                              : "reference box walk (GPW_CONTRACT_CUBE=0)")
                 <<";  task list: "<<itsBoxTasks.size()<<" shell pairs, "<<nTasks<<" (pair, offset) tasks, "
                 <<(double(nTasks*sizeof(BoxTask))/1048576.0)<<" MB"<<std::endl;
        if (qchem::report::Depth() > 0)
        {
            qchem::report::json s;
            s["kernel"]     = UseContractCube() ? "contract" : "walk";
            s["shellPairs"] = (long)itsBoxTasks.size();
            s["tasks"]      = (long)nTasks;
            s["bytes"]      = (long)(nTasks*sizeof(BoxTask));
            s["eps"]        = kDensityEps();
            qchem::report::EmitAt("grids", "boxTasks", s);
        }
        return itsBoxTasks;
    }
    //! The shell-pair visiting order for the THREADED loops: LONGEST FIRST, so a \c schedule(dynamic) loop
    //! cannot end on its biggest chunk (plan 3c reason 2).  The task COUNT is the proxy -- boxes within one
    //! shell pair differ only by the offset's prefactor, so the count tracks the work.
    //! \note The SERIAL loops keep the natural order: reordering changes the cross-pair grid reduction, and
    //! serial is the bit-anchor path.  A threaded run already reorders that reduction by construction.
    const std::vector<size_t>& BoxTaskOrder(const UnitCell& A) const
    {
        if (!itsBoxTaskOrder.empty()) return itsBoxTaskOrder;
        const std::vector<ShellPairTasks>& tl=BoxTasks(A);
        itsBoxTaskOrder.resize(tl.size());
        for (size_t i=0;i<tl.size();i++) itsBoxTaskOrder[i]=i;
        std::stable_sort(itsBoxTaskOrder.begin(),itsBoxTaskOrder.end(),
                         [&tl](size_t p, size_t q){return tl[p].tasks.size()>tl[q].tasks.size();});
        return itsBoxTaskOrder;
    }
    // --- T3 route (b) STREAM FOLD (doc/SymmetryUpgradePlan.md §6b) -----------------------------------------
    // Fold the (pair i<=j, offset n) collocation terms under the imposed crystal ops {W|tau}.  Pure
    // geometry-fixed bookkeeping, derived once like the box task list: per accepted op, the basis map
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
    //! (default \f$\Gamma\f$).  Returns the number of ops actually used.  Const like the box task list --
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
        // keeps the existing fold.  A CHANGE invalidates NOTHING but the fold itself: the box task list is
        // the eps-superset of every (shell pair, offset) term and the fold is applied as a PER-ITERATION
        // filter on top of it, so the same list serves a folded and a free run alike (it did not, while the
        // cache existed -- a REDUCED-built cache could never serve a free-run replay, §6b/T3.2).
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
        itsFoldReported=false;                            // a fold CHANGE gets its own [fold] line
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

    //! Announce the T3 pair-term fold in the shared [fold] format -- ALWAYS, including when it is off,
    //! because a disarmed run looks identical to an armed one on the console otherwise.  Emitted once, from
    //! the first collocation of a run (it used to ride the stream-cache build, which no longer exists).
    void ReportPairFold() const
    {
        if (itsFoldReported) return;
        itsFoldReported=true;
        const StreamFold* sf=itsStreamFold.get();
        const size_t nn=size();
        size_t nPairs=0, nRepPairs=0;
        for (auto i:indices()) for (auto j:indices(i))
        {
            nPairs++;
            if (!sf || (!sf->pairs[i*nn+j].dead && sf->pairs[i*nn+j].rep==int(i*nn+j))) nRepPairs++;
        }
        qchem::report::EmitFold("collocation pair terms (T3)", sf ? sf->maps.size() : 0,
                                nPairs, nRepPairs,
                                sf ? std::string() : std::string("free/multi-k run, or GPW_STREAM_FOLD=0"));
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
        // evaluates every (shell pair, offset) task live -- there is no value cache any more (plan step 7).
        // The TASK LIST build nested inside charges its own setup bucket, so this one reads as the pure
        // per-iteration scatter.
        qchem::report::Timed timer("scf: collocate density (pair scatter)");
        ReportPairFold();
        std::vector<rvec_t> rho(K);
        for (size_t l=0; l<K; l++) rho[l]=rvec_t(size_t(N_L[l].x)*N_L[l].y*N_L[l].z, 0.0);
        // Hermitian fold: the (j,i,-R) term contributes the SAME (idx, value, weight) as (i,j,R) -- the product
        // field is the lattice-translated twin and D_ji e^{-ik(-R)} = conj(D_ij e^{-ikR}) -- so loop j>=i and
        // double the off-diagonal weight (a 2x saving over the full n^2 loop, exact).
        const size_t nn=size();
        const StreamFold* sf=itsStreamFold.get();               // T3 route (b): reduced scatter (§6b), null = full
        // T3.5: the reduced scatter reads the ORBIT-PROJECTED D at each representative, never the
        // representative's own element -- see FoldProjectedD.  O(n^2) per call beside an O(pairs x points)
        // scatter, and it is what makes the folded density equal the unfolded imposed one for ANY iterate.
        const std::vector<dcmplx> Dsym=FoldProjectedD(D);
        const std::vector<ShellPairTasks>& tl=BoxTasks(A);
        // THE SHELL-BLOCKED SCATTER -- now the ONLY route.  Every component pair of ONE shell pair shares
        // its box walk: offsets, level, geometry, ellipsoid pre-screen and wrap are all shell properties.
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
            // The component pairs of this shell pair with a live D and (under a fold) a representative edge.
            std::vector<size_t> want;                                   // encoded i*nn+j
            for (size_t i=si.begin;i<si.end;i++)
                for (size_t j=std::max(i,sj.begin);j<sj.end;j++)
                {
                    if (sf) { const PairEdge& pe=sf->pairs[i*nn+j];
                              if (pe.dead || pe.rep!=int(i*nn+j)) continue; }
                    if ((sf ? Dsym[i*nn+j] : dcmplx(D(i,j)))==dcmplx(0.0)) continue;
                    want.push_back(i*nn+j);
                }
            if (want.empty()) return;
            // The level is a SHELL-pair property (PairLevel reads the shared radial's MaxExponent), asserted
            // per component below in a debug build.
            const size_t L=PairLevel(si.begin,sj.begin,ecut_L,0.0,0.0,relFieldSharp);
#ifndef NDEBUG
            for (size_t key : want)
                assert(L==PairLevel(key/nn,key%nn,ecut_L,0.0,0.0,relFieldSharp)
                       && "the pair->level assignment must be a shell-pair property");
#endif
            rvec_t& r=dst[L];
            std::vector<size_t> here;                                   // this offset's live pairs
            std::vector<double> cwHere, epsHere;
            // The component-LOCAL slots of each live pair in the shell-pair factor arrays.  Decoding them
            // from the packed key inside the point loop (key/nn, key%nn) put a HARDWARE INTEGER DIVIDE on
            // the innermost path -- 13.9% of the uncached kernel, measured 2026-08-26.  They are properties
            // of the (pair, offset) pre-filter, not of the point, so they are resolved once here.
            std::vector<size_t> iaHere, jbHere;
            for (const BoxTask& task : tl[sp].tasks)                     // the hoisted geometry (see BoxTasks)
            {
                const ivec3_t& n=task.n; const rvec3_t& Roff=task.Roff;
                here.clear(); cwHere.clear(); epsHere.clear(); iaHere.clear(); jbHere.clear();
                double epsUnion=0.0;
                const double pf=task.pf;                                 // screen (1), shared by the shell pair
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
                    const double e=UseDAwareScreen() ? std::max(kDensityEps(), kDensityEps()/std::fabs(c))
                                                     : kDensityEps();
                    if (pf>=-std::log(e)) continue;                     // whole box below THIS pair's eps
                    here.push_back(key); cwHere.push_back(c*wm); epsHere.push_back(e);
                    iaHere.push_back(i-si.begin); jbHere.push_back(j-sj.begin);
                    epsUnion = epsUnion==0.0 ? e : std::min(epsUnion,e);
                }
                if (here.empty()) continue;
                // THE SEPARABLE-CONTRACTION ROUTE (GPW_CONTRACT_CUBE, doc/CollocationRewritePlan.md step 5).
                // The whole shell pair collapses into ONE Gaussian x ONE polynomial with these weights
                // summed in, so the per-component-pair loop below disappears entirely.  Note the
                // per-component |v|<epsHere screen has NO analogue here and needs none: the box is already
                // sized by epsUnion, so that test never shrank the walk -- it only declined to accumulate a
                // value it had already computed (step 3).  Not accumulating is the only thing it bought,
                // and the contracted cube is formed regardless, so dropping it costs nothing and removes a
                // truncation.  The kernel DECLINES contracted shells and lp>=kMaxPoly; those fall through.
                bool didContract=false;
                if (UseContractCube())
                {
                    const size_t nIs=si.end-si.begin, nJs=sj.end-sj.begin;
                    double wd[kMaxShell*kMaxShell];
                    for (size_t t=0;t<nIs*nJs;t++) wd[t]=0.0;
                    for (size_t k=0;k<here.size();k++) wd[iaHere[k]*nJs+jbHere[k]]=cwHere[k];
                    const BoxGeom bg=MakeBoxGeom(si.begin,sj.begin,Roff,A,N_L[L],epsUnion);
                    if (bg.live)
                    {
                        const PairPoly q=MakePairPoly(si.begin,nIs,sj.begin,nJs,Roff,wd);
                        if (q.live) { ContractCube(bg,q,N_L[L],r); didContract=true; }
                    }
                }
                if (didContract) continue;
                // ⛔ NO PER-COMPONENT |v| SCREEN HERE ANY MORE (2026-08-27, plan step 7).  It used to read
                // `if (fabs(v)<epsHere[k]) continue;`, and dropping it is what makes the walk and the
                // contracted cube truncate IDENTICALLY -- same box, same terms -- which is what keeps
                // collocate/integrate an EXACT ADJOINT on either route.  Two facts settle it: the screen
                // never shrank the walk (the box is already sized by epsUnion, so it only declined to
                // accumulate a value it had ALREADY computed -- plan step 3), and while the value cache
                // existed both directions replayed ONE frozen stream, so the asymmetry it creates was
                // invisible.  Deleting the cache exposed it at once: GPW.AnalyticIntegrateBackAdjoint and
                // the two XC finite-difference gates all sat ~1.2e-8 relative off on the walk, scaling
                // linearly with GPW_DENSITY_EPS (green again at 1e-12).  epsHere survives only as the
                // input to epsUnion, which is where a component pair's own tolerance belongs.
                ForShellPairBox(si.begin,si.end-si.begin,sj.begin,sj.end-sj.begin,Roff,A,N_L[L],
                                [&](size_t idx, const double* fI, const double* fJ)
                {
                    for (size_t k=0;k<here.size();k++) r[idx]+=cwHere[k]*fI[iaHere[k]]*fJ[jbHere[k]];
                }, epsUnion);
            }
        };
#ifdef QCHEM_OPENMP
        if (const int nthreads=PairThreads(); nthreads>1)
        {
            // OpenMP over the SHELL pairs.  Each thread owns a PRIVATE level-density accumulator `mine`
            // (declared inside the parallel region -> genuinely thread-local, so no omp_get_thread_num /
            // <omp.h> needed), scatters its share into it, then folds it into the shared rho under a
            // critical section (the "per-thread accumulators + reduce").  The reduction reorders the
            // cross-pair grid sums, so a threaded run drifts a few ULPs from serial -- accepted; the Si
            // bit-anchors always run serial (GPW_OMP_THREADS unset -> nthreads==1 -> the branch below).
            // LONGEST FIRST (plan 3c reason 2): the task list already knows each shell pair's size, so the
            // dynamic loop cannot end on its biggest chunk.  Serial keeps the natural order (see BoxTaskOrder).
            const std::vector<size_t>& order=BoxTaskOrder(A);
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
                for (size_t q=0;q<order.size();q++)
                {
                    try { scatterShell(order[q],mine); }
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
        for (size_t sp=0;sp<sprs.size();sp++) scatterShell(sp,rho);       // serial (byte-identical default)
        return rho;
    }
    // --- Phase-independent integrate-back memo ---------------------------------------------------------------
    // h_ij(k) = w_l Sum_n e^{+ik.R_n} B_ij(n), with B_ij(n) = Sum_box chi_i chi_j^n V.  The expensive per-offset
    // reductions B are k-INDEPENDENT (the Bloch phase enters only the final contraction), and the SAME field V
    // is integrated repeatedly: the static local PP by EVERY k-block (once per SCF each), and the per-iteration
    // KS fields once per k-block within an iteration.  So memoize {B_ij(n)} keyed on the EXACT (ladder shape,
    // scale, V_L): the first caller pays the sweep, every other k-block only contracts phases (2026-07-15
    // multi-k profile: the per-k MakeLocalPP sweep alone was ~10% of the anchor).  Like the box task list, the
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
    // --- INTEGRATE-BACK CALL CENSUS (GPW_INTEGRATE_CENSUS=1) ----------------------------------------------
    // An INSTRUMENT, not a feature: off by default, costs one hash of the field per call when on.  It exists
    // because the collocation side's call count turned out to hide a 3.7x (the depth-1 CollocMemo, 2026-08-28)
    // and the same question is open here -- the per-iteration KS path passes `screenD`, and screened calls
    // bypass IntegrateMemo BY DESIGN (the active set moves with D while the memo is keyed on V alone), so
    // nothing has ever counted how many of those fields are actually distinct.
    //
    // THREE OUTCOMES, and they call for different fixes:
    //   NEW                      -- genuinely distinct work; the count is the physics.
    //   REPEAT (V and screen)    -- straight redundancy; a (V, screen)-keyed memo would replay it exactly.
    //   REPEAT (V only)          -- the SAME field re-integrated under a WIDENED screen.  Not replayable as
    //                              is, but it says the screen is what defeats the memo, which is a different
    //                              and more interesting problem than "the caller asks twice".
    static bool IntegrateCensus()
    { static const bool b=[]{const char* s=std::getenv("GPW_INTEGRATE_CENSUS"); return s && std::atoi(s)!=0;}(); return b; }
    //! FNV-1a over the raw bit patterns.  A diagnostic: a collision mislabels one line of a report and
    //! changes no result, which is why this is not the exact compare IntegrateMemo::SameField uses.
    static size_t BitHash(size_t h, double v)
    {
        unsigned long long b=0; std::memcpy(&b,&v,sizeof(b));
        for (int k=0;k<8;k++) { h^=size_t((b>>(8*k))&0xff); h*=1099511628211ull; }
        return h;
    }
    static size_t FieldHash(const std::vector<rvec_t>& V_L)
    {
        size_t h=1469598103934665603ull;
        for (const rvec_t& v : V_L) { h=BitHash(h,double(v.size())); for (double x : v) h=BitHash(h,x); }
        return h;
    }
    size_t ScreenHash(const chmat_t* screenD) const
    {
        if (!screenD) return 0;
        size_t h=1469598103934665603ull;
        for (size_t i=0;i<screenD->rows();i++)
            for (size_t j=i;j<screenD->columns();j++) h=BitHash(h,std::abs(dcmplx((*screenD)(i,j))));
        return h;
    }
    //! Classify this call against the recent history and charge the matching ledger bucket (the ledger's
    //! [xN] IS the census).  Keeps the last 32 (field, screen) hashes -- more than an SCF iteration's worth.
    void CensusIntegrateCall(const std::vector<rvec_t>& V_L, const chmat_t* screenD) const
    {
        const size_t fh=FieldHash(V_L), sh=ScreenHash(screenD);
        bool sameBoth=false, sameField=false;
        for (const auto& [f,c] : itsFieldHistory)
            if (f==fh) { sameField=true; if (c==sh) { sameBoth=true; break; } }
        qchem::report::Timed t(sameBoth  ? "census: h field REPEAT (V and screen seen)"
                             : sameField ? "census: h field REPEAT (V only -- screen widened)"
                                         : "census: h field NEW");
        itsFieldHistory.emplace_back(fh,sh);
        if (itsFieldHistory.size()>32) itsFieldHistory.erase(itsFieldHistory.begin());
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
        if (absRelCutoff==0.0 && !pairLevels)
        {
            timer.emplace("scf: integrate-back (pair gather)");
            if (IntegrateCensus()) CensusIntegrateCall(V_L, screenD);
        }
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
        // Miss: evaluate every (shell pair, offset) task live, recording B for the next caller.  There is
        // no value cache any more (plan step 7) -- the TASK LIST supplies the geometry and the contraction
        // kernel supplies the arithmetic, in both directions.
        const std::vector<ShellPairTasks>& tl=BoxTasks(A);
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
        // The SHELL-BLOCKED gather -- the exact mirror of `scatterShell`, and now the ONLY route.  It serves
        // both the per-iteration density path and the STATIC SHARP-FIELD sweeps (MakeLocalPP's kappa rule,
        // the explicit-pairLevels V_loc ball): every pair of a shell pair shares one box walk.  The D-aware
        // box resolves as the UNION of the component pairs' tolerances, exactly as on the density side, and
        // NEITHER side applies a per-component |v| screen on top of it -- which is what keeps collocate and
        // integrate on the SAME active set, hence exact adjoints (see the note in `scatterShell`).
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
                    if (sf) { const PairEdge& pe=sf->pairs[i*nn+j];
                              if (pe.dead || pe.rep!=int(i*nn+j)) continue; }
                    want.push_back(i*nn+j);
                }
            if (want.empty()) return;
            // The level is a shell-pair property in both regimes: the explicit pairLevels assignment
            // (StaticFieldPairLevels keys on MaxExponent(i)+MaxExponent(j)) and PairLevel itself both read
            // the shared radial.
            const size_t i0=want.front()/nn, j0=want.front()%nn;
            const size_t l = pairLevels ? (*pairLevels)[i0*nn+j0]
                           :              PairLevel(i0,j0,ecut_L,absRelCutoff,fieldSharpness,relFieldSharp);
            const rvec_t& V=V_L[l];
            const double  w=A.GetCellVolume()/double(V.size());
#ifndef NDEBUG
            for (size_t key : want)                              // every component pair must agree on it
                assert(l==(pairLevels ? (*pairLevels)[key]
                          :             PairLevel(key/nn,key%nn,ecut_L,absRelCutoff,fieldSharpness,relFieldSharp))
                       && "the pair->level assignment must be a shell-pair property");
#endif
            if (memo) for (size_t key : want) memo->B[key].level=l;
            std::vector<dcmplx> s(want.size(), dcmplx(0.0));
            std::vector<size_t> here;
            std::vector<double> wmHere, epsHere, bHere;
            std::vector<size_t> iaHere, jbHere;      // the component-LOCAL slots -- see scatterShell
            for (const BoxTask& task : tl[sp].tasks)                     // the hoisted geometry (see BoxTasks)
            {
                const ivec3_t& n=task.n; const rvec3_t& Roff=task.Roff;
                here.clear(); wmHere.clear(); epsHere.clear(); iaHere.clear(); jbHere.clear();
                double epsUnion=0.0;
                const double pf=task.pf;                                 // screen (1), shared by the shell pair
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
                    if (screenD && UseDAwareScreen())
                    {
                        // ⛔ THE SCREEN IS THE MAGNITUDE |D_ij|, NEVER Re[D conj(phase)] -- the SAME defect
                        // this direction's per-pair path was fixed for, and this copy of it survived because
                        // until 2026-08-27 the value cache sent every in-budget pair down the other path
                        // (doc/Benchmark.md footnote 1, the shifted-MP defect; deleting the cache in plan
                        // step 7 routed all of them here and it fired at once: Si 2x2x2 shifted MP came out
                        // 4.1 Ha high and NON-AUFBAU).  A real part is not a magnitude: this direction's
                        // term is phase(n)*b, whose size is |b| however the phase is oriented, so a term
                        // with Re[D conj(phase)] = 0 still adds a purely IMAGINARY, entirely non-negligible
                        // amount to H.  At a QUARTER-INTEGER k that is systematic, not accidental --
                        // e^{2 pi i k n} = i^n is purely imaginary for every ODD offset, and screenD is a
                        // matrix of MAGNITUDES (real, positive), so the old test discarded EVERY odd-offset
                        // term and H came out exactly real.  |D_ij| = |D_ij conj(phase)| is strictly MORE
                        // conservative, so it can only keep terms the old test dropped.
                        const double Dmag = ((i==j)?1.0:2.0)
                                          * (Dscreen.empty() ? std::abs(dcmplx((*screenD)(i,j)))
                                                             : Dscreen[key]);
                        if (Dmag==0.0) continue;
                        e=std::max(kDensityEps(), kDensityEps()/Dmag);
                    }
                    if (pf>=-std::log(e)) continue;                     // whole box below THIS pair's eps
                    here.push_back(k); wmHere.push_back(wm); epsHere.push_back(e);
                    iaHere.push_back(i-si.begin); jbHere.push_back(j-sj.begin);
                    epsUnion = epsUnion==0.0 ? e : std::min(epsUnion,e);
                }
                if (here.empty()) continue;
                bHere.assign(here.size(),0.0);
                // THE GATHER, the exact transpose of the scatter (doc/CollocationRewritePlan.md step 6).
                // Same tables, same chord, same coefficients -- so with GPW_CONTRACT_CUBE on BOTH directions
                // the collocate/integrate pair stays an exact adjoint.  Switching only one breaks it (that
                // is what moved Si's Etot by 2.4e-7 Ha in step 5a), so the two are deliberately wired to the
                // SAME flag.
                bool didContract=false;
                if (UseContractCube())
                {
                    const size_t nIs=si.end-si.begin, nJs=sj.end-sj.begin;
                    const BoxGeom bg=MakeBoxGeom(si.begin,sj.begin,Roff,A,N_L[l],epsUnion);
                    if (bg.live)
                    {
                        const PairPoly pq=MakePairPoly(si.begin,nIs,sj.begin,nJs,Roff,nullptr);
                        if (pq.live)
                        {
                            GridPoly W;
                            GatherCube(bg,pq,N_L[l],V,W);
                            double bAll[kMaxShell*kMaxShell];
                            MomentsToPairs(bg,pq,W,si.begin,nIs,sj.begin,nJs,bAll);
                            for (size_t qq=0;qq<here.size();qq++) bHere[qq]=bAll[iaHere[qq]*nJs+jbHere[qq]];
                            didContract=true;
                        }
                    }
                }
                if (!didContract)                    // no per-component |v| screen -- see scatterShell
                ForShellPairBox(si.begin,si.end-si.begin,sj.begin,sj.end-sj.begin,Roff,A,N_L[l],
                                [&](size_t idx, const double* fI, const double* fJ)
                {
                    const double Vv=V[idx];
                    for (size_t q=0;q<here.size();q++) bHere[q]+=fI[iaHere[q]]*fJ[jbHere[q]]*Vv;
                }, epsUnion);
                for (size_t q=0;q<here.size();q++)
                {
                    if (memo) memo->B[want[here[q]]].nb.emplace_back(n,bHere[q]);
                    s[here[q]]+=phase(n)*(wmHere[q]*bHere[q]);
                }
            }
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
            const std::vector<size_t>& order=BoxTaskOrder(A);   // longest first (see BoxTaskOrder)
            std::exception_ptr firstEx;   // throw containment -- see CollocateDensity
            #pragma omp parallel for schedule(dynamic) num_threads(nthreads)
            for (size_t q=0;q<order.size();q++)
            {
                try { integrateShell(order[q]); }
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
        for (size_t sp=0;sp<sprs.size();sp++) integrateShell(sp);           // serial (byte-identical default)
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
    mutable std::vector<std::pair<size_t,size_t>> itsFieldHistory;  //!< GPW_INTEGRATE_CENSUS (field, screen) hashes
    mutable std::vector<ShellPairTasks> itsBoxTasks;      //!< the (shell pair, offset) geometry, derived once
    mutable std::vector<size_t>        itsBoxTaskOrder;  //!< its longest-first permutation (threaded loops)
    mutable bool                       itsFoldReported=false;  //!< the [fold] line is emitted once per run
    mutable std::vector<IntegrateMemo> itsIntegrateMemos;  //!< phase-independent B_ij(n) per exact field (LRU)
    mutable std::unique_ptr<StreamFold> itsStreamFold;     //!< T3 route (b) (pair, R) orbit fold (null = free run)
};

static_assert(is1E_Evaluator <NR_Evaluator>, "NR_Evaluator must satisfy is1E_Evaluator");
static_assert(isDFT_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isDFT_Evaluator (ThreeC kernel)");
static_assert( isHF_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isHF_Evaluator (FourC kernel)");

} //namespace

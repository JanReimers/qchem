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
#include <vector>
#include <cmath>     // std::sqrt/std::log (the lattice-sum magnitude-screening reach radius)
#include <algorithm> // std::min (alpha_min per radial)
export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;
import qchem.BasisSet.Molecule.Evaluators;                             // Evaluator + concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;      // PGData
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;  // GaussianRF named kernels
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;// Polarization
import qchem.Structure;
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

    // Molecule::LatticeSum1E: the finest exponent over every component's radial primitives -- the density
    // grid's resolution driver (see LatticeSum1E::MaxExponent).  Walks the owned radials; no primitive escapes.
    double MaxExponent() const
    {
        double amax=0.0;
        for (auto i:indices()) for (double e : radials[i]->GetExponents()) if (e>amax) amax=e;
        return amax;
    }

    // --- GPW collocation adjoint (integrate-back), PATCHED (Molecule::LatticeSum1E::MakePotentialMatrix) -----
    // Per-orbital COMPRESSED COLUMNS of the Bloch orbital chi_i^k on the density grid: itsSupportIdx[i] lists
    // the grid-point indices where |chi_i^k| exceeds kScreenEps of its peak, and itsSupportVal[i] the matching
    // chi_i^k values (complex; real at Gamma).  A pair (i,j) then contracts M_ij = w Sum_p conj(chi_i^k) V
    // chi_j^k over the INTERSECTION of the two index lists (both ascending -> a merge) -- the sub-eps grid points
    // of EITHER orbital contribute < eps to the product and are dropped, so this is bit-consistent with the dense
    // Phi^H (V.*Phi) contraction to the screening tolerance.  Geometry-fixed, so built once and reused across SCF
    // iterations (V changes, the supports do not).  The Gaussian primitives stay encapsulated: only chi values +
    // grid-point indices are stored.  (Single-grid scaffold for the multi-grid rewrite -- see LatticeSum1E.)
    // Build one orbital's COMPRESSED COLUMN on a grid: the grid-point indices where the Bloch orbital
    // chi_i^k(r_p)=Sum_R e^{ik.R} chi_i(r_p-R) exceeds kScreenEps of its peak, plus the matching values.  The
    // image shift R is applied as evaluating the home component at q=r_p-R (the Gaussian depends on |r-centre|).
    void BuildColumn(const rvec3vec_t& gridPts, const std::vector<rvec3_t>& Rs, const cvec_t& phases,
                     size_t i, std::vector<int>& idx, std::vector<dcmplx>& val) const
    {
        idx.clear(); val.clear();
        const GaussianRF& rf=*radials[i];
        const rvec3_t ci=rf.GetCenter();
        const size_t Np=gridPts.size();
        std::vector<dcmplx> col(Np);
        double peak=0.0;
        for (size_t p=0;p<Np;p++)
        {
            dcmplx v(0.0);
            for (size_t r=0;r<Rs.size();r++)
            {
                const rvec3_t q=gridPts[p]-Rs[r];
                v += phases[r] * (ns[i]*pols[i](q-ci)*rf(q));
            }
            col[p]=v; peak=std::max(peak, std::abs(v));
        }
        const double thr=kScreenEps*peak;
        for (size_t p=0;p<Np;p++) if (std::abs(col[p])>thr) { idx.push_back(int(p)); val.push_back(col[p]); }
    }
    void EnsureSupports(const rvec3vec_t& gridPts, const std::vector<rvec3_t>& Rs, const cvec_t& phases) const
    {
        assert(phases.size()==Rs.size());
        if (itsSupportsBuilt && itsSupportNpts==gridPts.size()) return;   // built once per grid
        const size_t n=size();
        itsSupportIdx.assign(n,{}); itsSupportVal.assign(n,{});
        for (auto i:indices()) BuildColumn(gridPts, Rs, phases, i, itsSupportIdx[i], itsSupportVal[i]);
        itsSupportNpts=gridPts.size(); itsSupportsBuilt=true;
    }
    // Contract one pair over the OVERLAP of its two compressed columns (both ascending -> a merge):
    // s = Sum_p conj(chi_i^k) V[p] chi_j^k over p in supp(i) INTERSECT supp(j).  The sub-eps points of either
    // orbital contribute < eps to the product and are dropped -> bit-consistent with the dense contraction.
    static dcmplx ContractPair(const std::vector<int>& Ii, const std::vector<dcmplx>& Vi,
                               const std::vector<int>& Ij, const std::vector<dcmplx>& Vj, const rvec_t& V)
    {
        dcmplx s(0.0); size_t a=0,b=0;
        while (a<Ii.size() && b<Ij.size())
        {
            if      (Ii[a]<Ij[b]) ++a;
            else if (Ii[a]>Ij[b]) ++b;
            else { s += std::conj(Vi[a])*V[Ii[a]]*Vj[b]; ++a; ++b; }   // conj(chi_i^k) V chi_j^k (bra i conj)
        }
        return s;
    }
    chmat_t MakePotentialMatrix(const rvec3vec_t& gridPts, const std::vector<rvec3_t>& Rs,
                                const cvec_t& phases, const rvec_t& V, double w) const
    {
        EnsureSupports(gridPts, Rs, phases);
        const size_t n=size();
        chmat_t M(n);
        for (auto i:indices()) for (auto j:indices(i))   // j>=i (Hermitian upper triangle)
        {
            dcmplx s=w*ContractPair(itsSupportIdx[i],itsSupportVal[i], itsSupportIdx[j],itsSupportVal[j], V);
            M(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;   // Hermitian diagonal real; (j,i) mirrored as conj
        }
        return M;
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

    // --- MULTI-GRID patched integrate-back (Molecule::LatticeSum1E::MakePotentialMatrixMG, GPWPlan.md S0 Inc 2) -
    // Each orbital PAIR is contracted on the COARSEST grid LEVEL that still resolves its product exponent
    // alpha_i+alpha_j, so diffuse pairs live on small coarse grids instead of the fine grid dictated by the
    // TIGHTEST primitive: O(n^2 Npts_fine) -> O(sum_L pairs(L) Npts(L)).  The levels arrive finest-first --
    // gridPts_L[L] the level's grid points, V_L[L] the potential restricted to that level (low-passed by the
    // level's own {G}), w_L[L]=Omega/Npts(L) its quadrature weight, ecut_L[L] its cutoff (DESCENDING).  The
    // pair->level assignment uses the internal exponents (primitives stay encapsulated).  APPROXIMATE vs the
    // single fine grid (a coarse-level pair carries that level's adequate-by-construction grid error, CP2K's
    // REL_CUTOFF multigrid); converges to the fine-grid result as the ladder refines.  Per-(level,orbital)
    // compressed columns are built LAZILY (only where an orbital participates) and cached across SCF iterations.
    chmat_t MakePotentialMatrixMG(const std::vector<rvec3vec_t>& gridPts_L, const std::vector<double>& ecut_L,
                                  const std::vector<rvec3_t>& Rs, const cvec_t& phases,
                                  const std::vector<rvec_t>& V_L, const std::vector<double>& w_L) const
    {
        const size_t n=size(), K=gridPts_L.size();
        assert(K>0 && ecut_L.size()==K && V_L.size()==K && w_L.size()==K);
        EnsureMGCache(K, gridPts_L[0].size());
        const double amax=MaxExponent(), efine=ecut_L[0];
        std::vector<double> amaxi(n);                          // per-orbital tightest primitive (product-resolution)
        for (auto i:indices()) { double a=0.0; for (double e:radials[i]->GetExponents()) a=std::max(a,e); amaxi[i]=a; }
        chmat_t M(n);
        for (auto i:indices()) for (auto j:indices(i))
        {
            const double req=efine*(amaxi[i]+amaxi[j])/(2.0*amax);   // cutoff that resolves the pair product
            size_t L=0; for (size_t l=0;l<K;l++) if (ecut_L[l]>=req) L=l;  // coarsest (largest l) still >= req
            BuildMGColumn(L, i, gridPts_L[L], Rs, phases);
            BuildMGColumn(L, j, gridPts_L[L], Rs, phases);
            dcmplx s=w_L[L]*ContractPair(itsMGIdx[L*n+i],itsMGVal[L*n+i], itsMGIdx[L*n+j],itsMGVal[L*n+j], V_L[L]);
            M(i,j) = (i==j) ? dcmplx(std::real(s),0.0) : s;
        }
        return M;
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
    // (Re)size the per-(level,orbital) MG cache to K levels x n orbitals, clearing the built flags when the
    // level count or the fine-grid size changes (a new geometry).  Columns are then filled lazily by BuildMGColumn.
    void EnsureMGCache(size_t K, size_t nptsFine) const
    {
        if (itsMGK==K && itsMGNptsFine==nptsFine && !itsMGBuilt.empty()) return;
        const size_t nn=K*size();
        itsMGIdx.assign(nn,{}); itsMGVal.assign(nn,{}); itsMGBuilt.assign(nn,0);
        itsMGK=K; itsMGNptsFine=nptsFine;
    }
    // Build (once) the compressed column of orbital i on level L's grid; cached for reuse across SCF iterations.
    void BuildMGColumn(size_t L, size_t i, const rvec3vec_t& gridPts,
                       const std::vector<rvec3_t>& Rs, const cvec_t& phases) const
    {
        const size_t key=L*size()+i;
        if (itsMGBuilt[key]) return;
        BuildColumn(gridPts, Rs, phases, i, itsMGIdx[key], itsMGVal[key]);
        itsMGBuilt[key]=1;
    }

    // Per-orbital compressed columns of chi_i^k on the GPW density grid (see MakePotentialMatrix); built once
    // (geometry-fixed) by EnsureSupports, reused every SCF iteration.  Empty unless the GPW integrate-back runs.
    mutable std::vector<std::vector<int>>    itsSupportIdx;//!< grid-point indices where |chi_i^k| > eps*peak
    mutable std::vector<std::vector<dcmplx>> itsSupportVal;//!< the matching chi_i^k values (complex; real at Gamma)
    mutable size_t                        itsSupportNpts=0;//!< grid size the supports were built for (cache guard)
    mutable bool                          itsSupportsBuilt=false;
    // MULTI-GRID cache: per (level L, orbital i) compressed column, indexed [L*n+i]; built lazily by BuildMGColumn.
    mutable std::vector<std::vector<int>>    itsMGIdx;
    mutable std::vector<std::vector<dcmplx>> itsMGVal;
    mutable std::vector<char>                itsMGBuilt;   //!< per-(L,i) built flag
    mutable size_t                           itsMGK=0;     //!< level count the cache is sized for
    mutable size_t                           itsMGNptsFine=0; //!< fine-grid size (geometry-change guard)
};

static_assert(is1E_Evaluator <NR_Evaluator>, "NR_Evaluator must satisfy is1E_Evaluator");
static_assert(isDFT_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isDFT_Evaluator (ThreeC kernel)");
static_assert( isHF_Evaluator<NR_Evaluator>, "NR_Evaluator must satisfy isHF_Evaluator (FourC kernel)");

} //namespace

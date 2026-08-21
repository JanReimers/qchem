// File: ExactIrrepCD.C  Exact implementation of the charged density.
module;
#include <algorithm>   // std::min (the threaded rho-sampling block partition)
#include <cassert>
#include <complex>
#include <exception>   // std::exception_ptr (throw containment across the threaded block loop)
#include <iostream>
#include <stdexcept>
#include <stdlib.h>
#include <type_traits>
#include <vector>
#include <atomic>
#include <map>
#include <string>

module qchem.ChargeDensity.Imp.IrrepCD;
import qchem.Symmetry;
import qchem.Blaze;
import qchem.Parallel;                  // WorkerThreads (GPW_OMP_THREADS -- the rho-sampling GEMM)
import qchem.BasisSet.Orbital_DFT_IBS;   // cast the basis UP to the G-space capability (dcmplx path)

namespace qchem::ChargeDensity
{

// The density-freshness serial source (Version()) now lives in qchem.ChargeDensity as a SINGLE shared
// counter (NextDensityVersion): every density kind must draw from the same global clock or their serials
// collide across kinds and a dynamic term reuses a stale cached matrix (see the module doc + the SAD seed
// types NumericCD/SeedCD, which also stamp from it).

typedef Vector3D<std::complex<double> > Vec3;

rvec3_t  GradientContraction(const vec_t<rvec3_t >&, const vec_t<double>&, const rsmat_t&);
// rvec3_t  GradientContraction(const Vector<Vec3 >&, const Vector<std::complex<double> >&, const smat_t<std::complex<double> >&);

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> IrrepCD_Core<T>::IrrepCD_Core()
    : itsVersion(NextDensityVersion())
{};

template <class T> IrrepCD_Core<T>::IrrepCD_Core(const DenSMat& D,const tobs_t<T>* theBasisSet,Irrep qns)
    : itsDensityMatrix(D)
    , itsBasisSet(theBasisSet)
    , itsSpin(qns.ms)
    , itsIrrep(qns)
    , itsVersion(NextDensityVersion())
{
    assert(itsBasisSet);
};

template <class T> bool IrrepCD_Core<T>::IsZero() const
{
    // return max(abs(itsDensityMatrix))==0.0;
    return blazem::isZero(itsDensityMatrix);
}


//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
template <> void IrrepCD<double>::AccumulateDirect(rsmat_t& Jii) const
{
    const rohfbs_t* bs=dynamic_cast<const rohfbs_t*>(itsBasisSet);
    assert(bs);
    if (!IsZero()) bs->AccumulateDirect(Jii,itsDensityMatrix,bs);   // diagonal: this block's basis on both sides
}

template <> void IrrepCD<double>::AccumulateExchange(rsmat_t& Kii) const
{
    const rohfbs_t* bs=dynamic_cast<const rohfbs_t*>(itsBasisSet);
    assert(bs);
    if (!IsZero()) bs->AccumulateExchange(Kii,itsDensityMatrix,bs);
}

// One canonical irrep pair (this=i, other=j, i<=j) scattered into both irreps' Coulomb blocks.  This is the
// SOLE entry point the composite loops over (i<=j, so the DIAGONAL i==j lands here too): when other IS this
// (the self-pair) there is no bra-ket partner, so it collapses to a single localized contraction --
// ScatterBoth on the diagonal would add J.D + J^T.D = 2 J.D (the block is bra-ket symmetric).  Off-diagonal:
// the partner density is reached by a same-class cast (the IrrepCD<->IrrepCD idiom used by MixIn /
// GetChangeFrom); both empty -> nothing to build; the basis fetches ONLY the canonical J(i,j) block, so
// J(j,i) is never materialized.
template <class Leaf> void IrrepCD_HFPair<Leaf>::AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const tHF_Pair_CD<double>& other) const
{
    if (&other==this) { self().AccumulateDirect(Ji); return; }  // diagonal self-pair: single-block contraction
    // V1.8: no cast to the partner's concrete type.  We hand IT our block and let it finish with its own --
    // the contraction belongs to the BASIS, and a partner's only job is to present its matrix to that
    // routine.  Asking "give me your density matrix" would have reintroduced the getter IrrepCD is
    // deliberately without (CLAUDE.md cites this very class).
    const rohfbs_t* bs_i=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_i && "IrrepCD::AccumulateDirectBoth: this block's basis has no HF contraction face");
    other.CompleteDirectPair(Ji,Jj,self().itsDensityMatrix,bs_i);
}

// The partner's half: our own block completes the caller's pair.  Ji/Jj and Di/bs_i belong to the CALLER
// (irrep i); this object is irrep j.  The zero pre-screen lives here because only now are both blocks known.
template <class Leaf> void IrrepCD_HFPair<Leaf>::CompleteDirectPair(rsmat_t& Ji, rsmat_t& Jj, const rsmat_t& Di,
                                                    const BasisSet::Orbital_HF_IBS<double>* bs_i) const
{
    if (blazem::isZero(Di) && self().IsZero()) return;          // both empty: nothing to build
    const rohfbs_t* bs_j=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_j && "IrrepCD::CompleteDirectPair: this block's basis has no HF contraction face");
    bs_i->AccumulateDirectBoth(Ji,Jj,Di,self().itsDensityMatrix,bs_j);
}

// --- V1.31 whole-system route ---------------------------------------------------------------------------
// The basis answers whether it has one; this block only has to fold its own density up to the AO space.  The
// SLICE back down needs no density, so the composite drives that through the basis face directly.
template <> const BasisSet::WholeSystemFock_IBS<double>* IrrepCD<double>::WholeSystemFock() const
{
    return dynamic_cast<const BasisSet::WholeSystemFock_IBS<double>*>(itsBasisSet);
}

template <> void IrrepCD<double>::AddAODensity(rsmat_t& Dao) const
{
    const BasisSet::WholeSystemFock_IBS<double>* ws=WholeSystemFock();
    assert(ws && "IrrepCD::AddAODensity: this block's basis has no whole-system Fock route");
    if (!IsZero()) ws->AddAODensity(Dao,itsDensityMatrix);
}

// Exchange counterpart (see AccumulateDirectBoth): diagonal -> single AccumulateExchange, else the canonical
// Exchange block is fetched once.
template <class Leaf> void IrrepCD_HFPair<Leaf>::AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const tHF_Pair_CD<double>& other) const
{
    if (&other==this) { self().AccumulateExchange(Ki); return; }  // diagonal self-pair: single-block contraction
    const rohfbs_t* bs_i=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_i && "IrrepCD::AccumulateExchangeBoth: this block's basis has no HF contraction face");
    other.CompleteExchangePair(Ki,Kj,self().itsDensityMatrix,bs_i);
}

template <class Leaf> void IrrepCD_HFPair<Leaf>::CompleteExchangePair(rsmat_t& Ki, rsmat_t& Kj, const rsmat_t& Di,
                                                      const BasisSet::Orbital_HF_IBS<double>* bs_i) const
{
    if (blazem::isZero(Di) && self().IsZero()) return;
    const rohfbs_t* bs_j=dynamic_cast<const rohfbs_t*>(self().itsBasisSet);
    assert(bs_j && "IrrepCD::CompleteExchangePair: this block's basis has no HF contraction face");
    bs_i->AccumulateExchangeBoth(Ki,Kj,Di,self().itsDensityMatrix,bs_j);
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
// AO density-fit projection <rho|c> = Sum_ab D_ab <ab|c>, the finite (double) path's ProjectedDensity_AO
// face.  The periodic (dcmplx) density is NOT a ProjectedDensity_AO (see ProjectedDensityBase), so this is
// never reached for dcmplx; the if-constexpr keeps the double-only 3-centre machinery out of that build.
template <class T> rvec_t IrrepCD<T>::GetRepulsion3C(const BasisSet::rFIT_CD_ABS* fbs) const
{
    if constexpr (std::is_same_v<T,double>)
    {
        if (this->IsZero()) return rvec_t(fbs->GetNumFunctions(),0.0);
        auto dftbs=dynamic_cast<const todftbs_t<T>*>(this->itsBasisSet);
        assert(dftbs);
        // Contract the density matrix against the basis's CACHED, D-free 3-centre projection tensor <ab|c>
        // HERE -- the DENSITY owns D, so the D-contraction is a density operation, not a basis one.  The
        // basis exposes only the tensor (Repulsion3C(c) -> Projector3, built once, keyed by BasisSetID);
        // D never crosses into qcBasisSet.  This is the real-space model for fixing MakeFourierDensity(D): the
        // {G} 3-centre integral is the delta <ij|Dm>, and rho-tilde = Sum_ij D_ij <ij|Dm> is the SAME contraction.
        const auto& R=dftbs->Repulsion3C(*fbs).dense;     // <ab|c> (dense realization: one smat per fit function c)
        rvec_t ret(fbs->GetNumFunctions());
        for (size_t i=0;i<R.size();++i)
            ret[i]=blazem::sum(this->itsDensityMatrix % R[i]);  // <rho|c_i> = Sum_ab D_ab <ab|c_i>
        return ret;
    }
    else
        return rvec_t();   // inert: a periodic density carries no AO projection
}


// Energy = trace(D V) = Sum_ij D_ij V_ji = Sum_ij D_ij trans(V)_ij = sum(D % trans(V)).  For real
// symmetric V this is identical to sum(D % V) (the old form); for complex HERMITIAN D,V the transpose
// matters -- sum(D % V) would silently give the wrong (still-real) value.
template <class T> double IrrepCD_Core<T>::DM_Contract(const tStatic_CC<T>* v) const
{
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(v->GetMatrix(itsBasisSet,itsSpin)));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

template <class T> double IrrepCD_Core<T>::DM_Contract(const tDynamic_CC<T>* v,const tDM_CD<T>* cd) const
{
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(v->GetEMatrix(itsBasisSet,itsSpin,cd)));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

// This irrep's contribution to a whole-system energy: sum(D % B_i^T) with B_i the block for this basis.
// Same contraction as DM_Contract above, but the Fock block comes from a caller-supplied map (the HF
// term's cached blocks) instead of a per-irrep GetMatrix call -- so no GetMatrix round-trip.
template <class T> double IrrepCD_Core<T>::DM_ContractBlocks(const std::map<std::string,hmat_t<T>>& blocks) const
{
    if (IsZero()) return 0.0;
    T ComplexE=blazem::sum(itsDensityMatrix % blazem::trans(blocks.at(itsBasisSet->BasisSetID())));
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

// This irrep's rho at the caller's points from its pre-evaluated basis table (the FIELD dual of
// DM_ContractBlocks above): rho_g = Re[Phi D Phi^dag]_gg = Re sum_j (Phi D)_gj conj(Phi_gj) -- one
// (npts x n)(n x n) GEMM + a cheap rowwise dot, matching operator()(r)'s trans(phi) D conj(phi) per row.
// No table for this basis -> self-evaluate pointwise (correct; the caller's first pass may not have
// built every block's table yet).
// GPW_DM_RANK=1: the LOW-RANK CENSUS of D (doc/OpenWork.md Step 3, the rho-GEMM plan).  rho_g =
// Phi_g^dag D Phi_g costs O(npts n^2) here; if D = L L^dag with rank r, rho_g = ||L^dag Phi_g||^2 costs
// O(npts n r) and the win is n/r.  D is a density matrix so r SHOULD be the occupied count (~13/spin on
// MnO against n=118 => 6-9x), but SMEARING puts a fractional tail on the occupations, so the NUMERICAL
// rank at a usable tolerance is the quantity that decides it.  Measure before building any of it.
//
// Reports BOTH questions at once:
//   * lambda_min < 0 answers "is D PSD?".  A pivoted Cholesky needs PSD, and PSD is a property of the
//     MIXERS, not a theorem: linear mixing with alpha in [0,1] is CONVEX and preserves it, an
//     EXTRAPOLATION (density-space Pulay, alpha>1, MP/cold smearing's negative occupations) does not.
//     Today only LinearMixer touches D and its alpha is clamped <= 1, so this should read >= 0 -- and if
//     it ever does not, that IS the canary, not a surprise.
//   * the rank at a ladder of tolerances is the achievable n/r.
template <class T> void ReportDMRank(const hmat_t<T>& D, const Irrep& ir)
{
    static const bool on=std::getenv("GPW_DM_RANK")!=nullptr;
    if (!on || D.rows()==0) return;
    const size_t n=D.rows();
    rvec_t w; mat_t<T> U;
    blazem::eigen(D, w, U);                       // ascending; O(n^3), nothing beside the npts n^2 GEMM
    double tr=0.0; for (size_t i=0;i<n;i++) tr+=w[i];
    const double wmax=w[n-1];
    // PSD test on a RELATIVE floor, not on w[0]<0: an eigensolver returns O(eps*wmax) negative values for a
    // numerically-PSD matrix, so a strict <0 test cries wolf on every run (it did, first cut: lambda_min
    // ~ -1.8e-15 against wmax 13.2, i.e. ~1 ulp).  Below the floor = roundoff; a REAL negative eigenvalue
    // from an extrapolated D would sit orders of magnitude above it.
    const double psdFloor=-1e-12*wmax;
    std::cout<<"[DM rank] irrep="<<ir<<" n="<<n<<" Tr(D)="<<tr        // NB Tr(D), not Tr(DS)=N_elec
             <<" lambda: min="<<w[0]<<" ("<<w[0]/wmax<<" rel) max="<<wmax
             <<(w[0]<psdFloor ? "  ** NOT PSD **" : "  PSD(to roundoff)");
    for (double rel : {1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14})
    {
        size_t rk=0; for (size_t i=0;i<n;i++) if (w[i]>rel*wmax) rk++;
        std::cout<<"  r("<<rel<<")="<<rk;
    }
    std::cout<<std::endl;
    // THE EIGEN FACTOR'S LOCALITY, for side-by-side with the Cholesky factor's ([DM factor] below).
    // The eigenvectors are the NATURAL ORBITALS: the minimal-rank factorisation, but expected DELOCALISED.
    // The literature claim for pivoted Cholesky is the opposite -- localized MOs straight from D -- so this
    // pair of lines is the direct test of it, and it matters because THE TWO CONSUMERS WANT DIFFERENT
    // FACTORS: the rho GEMM is O(npts n r) so it wants the SMALLEST r (eigen), while collocation wants
    // compact boxes so it wants LOCALITY (Cholesky).  IPR is scale-free, so U vs U*sqrt(lambda) is moot.
    {
        const double keep=1e-10*wmax;
        double iprSum=0.0, iprMin=1e300, iprMax=0.0; size_t cnt=0; double ladder[5]={0,0,0,0,0};
        for (size_t k=0;k<n;k++)
        {
            if (w[k]<=keep) continue;                 // ascending: these are the surviving modes
            double s2=0.0,s4=0.0,cmax=0.0;
            for (size_t i=0;i<n;i++)
            { const double a=std::norm(std::complex<double>(U(i,k))); s2+=a; s4+=a*a; cmax=std::max(cmax,a); }
            const double ipr=(s4>0.0)?s2*s2/s4:0.0;
            iprSum+=ipr; iprMin=std::min(iprMin,ipr); iprMax=std::max(iprMax,ipr); cnt++;
            for (size_t t=0;t<5;t++)
            {
                const double f=std::pow(10.0,-2.0*double(t+1));
                size_t c=0; for (size_t i=0;i<n;i++)
                    if (std::norm(std::complex<double>(U(i,k)))>f*cmax) c++;
                ladder[t]+=double(c);
            }
        }
        // DOES LOCALITY CORRELATE WITH lambda? (user, 2026-08-21)  With lambda kept SEPARATE,
        // rho = sum_m lambda_m |phi_m|^2 is a SUM OF PER-MODE DENSITIES, so each mode can be collocated on
        // its OWN ladder level and the level densities summed -- exactly as pairs are.  If the modes'
        // spatial character tracks lambda, the level assignment is then free: read it off the spectrum.
        // Printed per mode, descending in lambda, so the correlation (or its absence) is visible directly.
        {
            std::cout<<"[DM lambda-vs-IPR]";
            for (size_t k=n; k-- > 0; )
            {
                if (w[k]<=keep) continue;
                double s2=0.0,s4=0.0;
                for (size_t i=0;i<n;i++)
                { const double a=std::norm(std::complex<double>(U(i,k))); s2+=a; s4+=a*a; }
                std::cout<<"  ("<<w[k]<<","<<((s4>0.0)?s2*s2/s4:0.0)<<")";
            }
            std::cout<<std::endl;
        }
        if (cnt>0)
            std::cout<<"[DM eigen] modes(1e-10)="<<cnt<<"  IPR: min="<<iprMin<<" mean="<<iprSum/double(cnt)
                     <<" max="<<iprMax
                     <<"  decay (mean #fns > 10^-t*max) t=1:"<<ladder[0]/double(cnt)
                     <<" t=2:"<<ladder[1]/double(cnt)<<" t=3:"<<ladder[2]/double(cnt)
                     <<" t=4:"<<ladder[3]/double(cnt)<<" t=5:"<<ladder[4]/double(cnt)<<std::endl;
    }
}

//! \brief LOW-RANK FACTOR of a density matrix: \f$D=LL^\dagger\f$ with \a L (n x rank), by PIVOTED
//! Cholesky.  Returns false -- leaving the caller on the full-rank path -- if D is not PSD.
//!
//! WHY THIS PAYS (doc/OpenWork.md Step 3).  \f$\rho_g=\Phi_g^\dagger D\,\Phi_g\f$ costs O(npts n^2); with
//! \f$D=LL^\dagger\f$ it is \f$\lVert L^\dagger\Phi_g\rVert^2\f$, i.e. O(npts n r).  D is a density matrix,
//! so its rank is the OCCUPIED count: MEASURED at 14-17 against n=118 on the MnO benchmark (GPW_DM_RANK=1)
//! => 7-8.4x, and the rank is the same from tol 1e-6 to 1e-12, so the occupied block is cleanly separated
//! from the null space and the cut is not a tuning decision.
//!
//! PIVOTED CHOLESKY, NOT A TRIMMED EIGENDECOMPOSITION (user pin).  Both give an L; the eigen route trims
//! clustered near-null eigenvectors, which are ill-conditioned (free to rotate within the cluster) and
//! carry that noise into the factor.  In the user's experience that damages LATE GDM CONVERGENCE.  The
//! objection is strictly weaker HERE, since we only ever MULTIPLY the factor and never invert it -- but
//! Cholesky is also ~3x cheaper and greedy on the diagonal, so there is no reason to prefer the other.
//!
//! EXACT, NOT AN APPROXIMATION: the pivot floor is LAPACK's own default (tol<0 => ~n*eps*max diag), so the
//! discarded residual sits AT ROUNDOFF.  This is not an accuracy trade like the Becke eps.
template <class T> bool LowRankFactor(const hmat_t<T>& D, mat_t<T>& L, size_t& rank)
{
    const size_t n=D.rows();
    if (n==0) return false;
    mat_t<T> Um(D);                                   // pstrf needs a plain (non-adaptor) matrix
    std::vector<blazem::blas_int_t> piv(n);
    const size_t m=(size_t)blazem::pstrf(Um, 'U', piv.data(), -1.0);   // P^T D P = U^H U, LAPACK default floor
    if (m==0 || m>n) return false;

    // L = P U^H, i.e. L(piv[a],k) = conj(U(k,a)) -- U is m x n upper TRAPEZOIDAL, so k<=a only, and pstrf
    // leaves the strict-lower part undefined (never read here).  Then L L^H = P U^H U P^T = D.
    L=mat_t<T>(n, m, T(0));
    for (size_t a=0;a<n;++a)
        for (size_t k=0;k<m && k<=a;++k)
            L(piv[a],k)=blazem::conj(Um(k,a));

    // PSD GUARD.  LAPACK pstrf does NOT fail on an indefinite matrix -- it simply stops when the pivot goes
    // non-positive and reports the rank it reached, so a non-PSD D would silently yield a TRUNCATED (wrong)
    // rho rather than an error.  Catch it on conserved mass: Tr(L L^H) = ||L||_F^2 must reproduce Tr(D).
    // A negative eigenvalue leaves that unaccounted for.  Cheap (O(n m)) beside the O(npts n r) GEMM.
    double trD=0.0;  for (size_t i=0;i<n;++i) trD+=blazem::real(D(i,i));
    double trL=0.0;  for (size_t i=0;i<n;++i) for (size_t k=0;k<m;++k) trL+=std::norm(std::complex<double>(L(i,k)));
    if (std::abs(trL-trD) > 1e-8*std::max(std::abs(trD),1.0)) return false;
    rank=m;

    // GPW_DM_RANK=1: ARE THE CHOLESKY ORBITALS LOCALIZED?  THE decisive unmeasured quantity for the
    // "factor D" idea (doc/OpenWork.md, its own section).  psi_m = sum_i L_im chi_i collocates over the
    // UNION of the boxes of the chi_i it actually touches -- so if each column of L is carried by a few
    // basis functions the orbitals have COMPACT boxes and the singles route can beat pair collocation; if
    // the columns are spread over all n, each psi_m costs the whole grid and that route is dead.
    // Measured centre-free by the INVERSE PARTICIPATION RATIO, IPR_m = (sum_i |L|^2)^2 / sum_i |L|^4 = the
    // EFFECTIVE number of basis functions carrying orbital m (n = fully delocalised, few = localised),
    // plus the count above a relative floor.  Pivoted Cholesky is greedy on the diagonal, so the
    // literature expectation ("Cholesky orbitals") is that these come out LOCALISED -- this measures it.
    static const bool on=std::getenv("GPW_DM_RANK")!=nullptr;
    if (on)
    {
        double iprMin=1e300, iprMax=0.0, iprSum=0.0; size_t sigMax=0, sigSum=0;
        double ladder[5]={0,0,0,0,0};
        for (size_t k=0;k<m;++k)
        {
            double s2=0.0, s4=0.0, cmax=0.0;
            for (size_t i=0;i<n;++i)
            {
                const double a=std::norm(std::complex<double>(L(i,k)));
                s2+=a; s4+=a*a; cmax=std::max(cmax,a);
            }
            const double ipr = (s4>0.0) ? s2*s2/s4 : 0.0;
            iprMin=std::min(iprMin,ipr); iprMax=std::max(iprMax,ipr); iprSum+=ipr;
            size_t sig=0; for (size_t i=0;i<n;++i)
                if (std::norm(std::complex<double>(L(i,k))) > 1e-6*cmax) sig++;   // 1e-3 in amplitude
            sigMax=std::max(sigMax,sig); sigSum+=sig;
            for (size_t t=0;t<5;++t)                      // the DECAY PROFILE: |L|>10^-(t+1) * max
            {
                const double f=std::pow(10.0,-2.0*double(t+1));   // in |L|^2
                size_t c=0; for (size_t i=0;i<n;++i)
                    if (std::norm(std::complex<double>(L(i,k)))>f*cmax) c++;
                ladder[t]+=double(c);
            }
        }
        std::cout<<"[DM factor] rank="<<m<<" of n="<<n
                 <<"  IPR(effective basis fns/orbital): min="<<iprMin<<" mean="<<iprSum/double(m)
                 <<" max="<<iprMax
                 <<"  |L|>1e-3*max count: mean="<<double(sigSum)/double(m)<<" max="<<sigMax
                 <<"\n[DM factor] coefficient DECAY (mean #fns with |L| > 10^-t * max), n="<<n
                 <<", ~"<<n/4<<" fns/atom:  t=1:"<<ladder[0]/double(m)<<"  t=2:"<<ladder[1]/double(m)
                 <<"  t=3:"<<ladder[2]/double(m)<<"  t=4:"<<ladder[3]/double(m)
                 <<"  t=5:"<<ladder[4]/double(m)<<std::endl;
    }
    return true;
}

template <class T> rvec_t IrrepCD_Core<T>::DM_RhoAtPoints(const rvec3vec_t& r, const std::map<Irrep,mat_t<T>>& Phi) const
{
    rvec_t ro(r.size(), 0.0);
    if (IsZero()) return ro;
    // Spatial key: Phi is a property of the basis block alone, so the two spin channels of one
    // spatial block share one table (Irrep interleaves spin into its ordering, hence Spin::None).
    auto it=Phi.find(Irrep(Spin::None,itsIrrep.sym));
    if (it==Phi.end())
    {
        for (size_t g=0; g<r.size(); g++) ro[g]=(*this)(r[g]);
        return ro;
    }
    const mat_t<T>& P=it->second;
    assert(P.rows()==r.size());
    assert(P.columns()==itsDensityMatrix.rows());
    // rho_g = Phi_g^dag D Phi_g, as (P D) row-dotted back into P: O(npts n^2), paid once per density
    // serial and THE largest per-iteration bucket of the atom-centred XC route (MnO Γ: 12 s/iteration
    // at 48k points x 122 functions).  Threaded by MESH-POINT BLOCK: rho_g depends on row g alone, so
    // each block is an independent GEMM + dot and the result is bit-identical at any thread count (a
    // partition of the OUTPUT, not of a reduction).  Opt-in via GPW_OMP_THREADS (qchem.Parallel).
    // D as a PLAIN matrix: an adaptor operand (HermitianMatrix) is not BLAS-dispatchable, so under
    // BLAZE_BLAS_MODE the product would silently fall back to blaze's own kernel -- 16x slower than
    // zgemm on this shape (measured).  One n x n copy per call is nothing beside the npts x n x n GEMM.
    ReportDMRank<T>(itsDensityMatrix, itsIrrep);   // hmat_t<T> is a conditional_t: T is non-deduced
    const mat_t<T> D(itsDensityMatrix);
    // LOW-RANK ROUTE.  D is a density matrix, so rank = the OCCUPIED count, not n: with D = L L^dag the
    // GEMM is (npts x n)(n x r) instead of (npts x n)(n x n), and rho_g = ||Psi_g||^2 replaces the row-dot.
    // MEASURED r=14-17 against n=118 on MnO => 7-8.4x ON THIS GEMM.  Exact to roundoff (LAPACK's own pivot
    // floor), so it is a cost change, not an accuracy trade -- verified by an A/B at identical energy.
    // ⚠ SCOPE, measured 2026-08-20: on the GPW *Kerker/Pulay* recipes this GEMM is nearly BYPASSED -- the
    // mixer hands XC a rho-tilde-backed density with NO D from iteration 1 on, so the hot path is the
    // matrix-free sampling in PWTerms, not this.  Measured 1.70 s here against 35.0 s there (6 iterations).
    // The 7-8x is therefore REAL but small on those recipes; it is worth having for the DM-backed routes
    // (molecular/atomic DFT, non-rho-tilde-mixed recipes) and grows with n.  Do NOT quote it as a
    // benchmark-row win.
    // GUARDED, not assumed: LowRankFactor returns false for a non-PSD D and we stay on the full path.  D is
    // PSD for any density assembled from orbitals with non-negative occupations -- which is every IrrepCD
    // today, because the only mixer that touches D is LINEAR with alpha in [0,1] (a CONVEX combination
    // preserves PSD) while Kerker/Pulay/Broyden extrapolate the FIELD, not D.  The fallback exists for the
    // paths that would break that: a density-space Pulay, alpha>1, MP/cold smearing's negative
    // occupations, or anything putting a DIFFERENCE of densities behind a DM face.
    // QCHEM_DM_LOWRANK=0 forces the full-rank path (the A/B valve; also the escape hatch if a future
    // density ever needs it).  2*rank < n is the pay-for-itself test: below that the thin GEMM plus the
    // O(n^3) factorisation is not obviously cheaper than the square one.
    static const bool lowRankOn=[]{ const char* e=std::getenv("QCHEM_DM_LOWRANK"); return !e || std::atoi(e)!=0; }();
    mat_t<T> L; size_t rank=0;
    const bool lowRank = lowRankOn && LowRankFactor<T>(itsDensityMatrix, L, rank) && 2*rank < P.columns();
    auto block=[&](size_t g0, size_t g1)
    {
        auto Pb = blazem::submatrix(P, g0, size_t(0), g1-g0, P.columns());
        if (lowRank)
        {
            mat_t<T> Psi = Pb*L;                          // (npts_b x rank) -- the thin GEMM
            for (size_t g=g0; g<g1; g++)
            {
                double acc=0.0;
                for (size_t m=0; m<Psi.columns(); m++) acc+=std::norm(std::complex<double>(Psi(g-g0,m)));
                ro[g]=acc;                                // rho_g = ||L^dag Phi_g||^2, real by construction
            }
            return;
        }
        mat_t<T> PD = Pb*D;
        for (size_t g=g0; g<g1; g++)
        {
            T acc{};
            for (size_t j=0; j<PD.columns(); j++)
                if constexpr (std::is_same_v<T,dcmplx>) acc+=PD(g-g0,j)*std::conj(P(g,j));
                else                                    acc+=PD(g-g0,j)*P(g,j);
            ro[g]=std::real(acc);
        }
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1 && r.size()>size_t(nthreads))
    {
        const size_t blk=(r.size()+size_t(nthreads)-1)/size_t(nthreads);
        std::exception_ptr firstEx;                   // throw containment (an escape = std::terminate)
        #pragma omp parallel for schedule(static,1) num_threads(nthreads)
        for (size_t b=0; b<(r.size()+blk-1)/blk; b++)
        {
            try { block(b*blk, std::min(r.size(), (b+1)*blk)); }
            catch (...)
            {
                #pragma omp critical (cd_rho_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return ro;
    }
#endif
    // Serial: multiply the WHOLE table, not a full-height view -- a submatrix operand costs blaze its
    // aligned/padded kernel choice, and this is the default (unthreaded) path every suite run takes.
    mat_t<T> PD = P*D;
    for (size_t g=0; g<r.size(); g++)
    {
        T acc{};
        for (size_t j=0; j<PD.columns(); j++)
            if constexpr (std::is_same_v<T,dcmplx>) acc+=PD(g,j)*std::conj(P(g,j));
            else                                    acc+=PD(g,j)*P(g,j);
        ro[g]=std::real(acc);
    }
    return ro;
}

template <class T> double IrrepCD_Core<T>::GetTotalCharge() const
{
    // N = integral rho = Sum_ij D_ij S_ji = Tr(D S) = sum(D % trans(S)) (% is the blaze Schur/direct product).
    // The trans matters for a genuinely COMPLEX Hermitian D,S (the complex-k Bloch case): sum(D % S) would be
    // Tr(D S^T) -- a different (still-real) value.  Same transpose fix as DM_Contract above; a no-op for real
    // symmetric S (trans(S)==S, byte-identical at Gamma / half-integer k).
    return std::real(blazem::sum(itsDensityMatrix%blazem::trans(itsBasisSet->Overlap())));
}


//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void IrrepCD_Core<T>::ReScale(double factor)
{
    // No UT coverage
    itsDensityMatrix*=factor;
    itsVersion=NextDensityVersion();   // a mutation is logically a new density
}

template <class T> void IrrepCD_Core<T>::MixIn(const tMixableDensity<T>& cd,double c)
{
    // Same-KIND cast to the CORE: both leaves (finite, periodic) mix identically -- D is D.  The
    // basis-ID assert below still pins that only matching blocks ever pair.
    const IrrepCD_Core<T>* eicd = dynamic_cast<const IrrepCD_Core<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    itsDensityMatrix = itsDensityMatrix*(1-c) + eicd->itsDensityMatrix*c;
    itsVersion=NextDensityVersion();   // a mutation is logically a new density
}

template <class T> double IrrepCD_Core<T>::GetChangeFrom(const tMixableDensity<T>& cd) const
{
    const IrrepCD_Core<T>* eicd = dynamic_cast<const IrrepCD_Core<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    return std::real(blazem::norm(itsDensityMatrix - eicd->itsDensityMatrix));
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double IrrepCD_Core<T>::operator()(const rvec3_t& r) const
{
    vec_t<T> phir=(*itsBasisSet)(r);
    return std::real(blazem::trans(phir)*itsDensityMatrix*blazem::conj(phir));
}

template <class T> rvec3_t IrrepCD<T>::Gradient(const rvec3_t& r) const
{
    // No UT coverage
    vec_t<T> phir=(*this->itsBasisSet)(r);
    vec_t<rvec3_t > gphir=this->itsBasisSet->Gradient(r);
    return GradientContraction(gphir,phir,this->itsDensityMatrix);
}

// rho-tilde from the density matrix, via the basis's G-space capability (plane-wave / dcmplx only).
// itsDensityMatrix already carries the BZ weight w_k (TOrbitals::GetChargeDensity scales it), so the
// composite's sum over blocks is the BZ average Sum_k w_k rho_k.
template <class Leaf> ΔG_Map IrrepCD_Fourier<Leaf>::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    // Contract D against the basis's D-free OVERLAP tensor (empty kernel) HERE -- the DENSITY owns D, so
    // this is where rho-tilde(dm) = (1/Omega) Sum_{G_i-G_j=dm} D_ij is formed.  The overlap-metric sibling
    // of GetRepulsion3C above; D never crosses into the basis.
    auto* fb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<typename Leaf::scalar_t,dcmplx>*>(self().itsBasisSet);
    assert(fb && "GetFourierDensity requires a periodic (G-space) Orbital_DFT_IBS basis");
    return Contract(fb->Overlap3C(c), self().itsDensityMatrix);
}

// Raw rho_DM(r) on c's raster (doc/GPWPlan 0.5(f2)): the basis's collocation-native forward, when it has
// one (GPW's Overlap3C carries applyRaw; a plane-wave basis leaves it empty -> we answer empty and the
// caller falls back to the ball route).  itsDensityMatrix carries the BZ weight, exactly as above.
template <class Leaf> rvec_t IrrepCD_Fourier<Leaf>::GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const
{
    auto* fb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<typename Leaf::scalar_t,dcmplx>*>(self().itsBasisSet);
    assert(fb && "GetRhoOnGrid requires a periodic (G-space) Orbital_DFT_IBS basis");
    const Projector3<dcmplx>& g=fb->Overlap3C(c);
    return g.applyRaw ? g.applyRaw(self().itsDensityMatrix) : rvec_t{};
}

// V_H(dm) = 4pi rho-tilde(dm)/|G|^2: contract D against the basis's D-free COULOMB tensor (Repulsion3C(c), the
// diagonal kernel baked in) -- the reciprocal mirror of the finite GetRepulsion3C(fbs) above.  D stays here.
template <class Leaf> ΔG_Map IrrepCD_Fourier<Leaf>::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    auto* fb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<typename Leaf::scalar_t,dcmplx>*>(self().itsBasisSet);
    assert(fb && "GetRepulsion3C requires a periodic (G-space) Orbital_DFT_IBS basis");
    return Contract(fb->Repulsion3C(c), self().itsDensityMatrix);
}


//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& IrrepCD_Core<T>::Write(std::ostream& os) const
{
    return os << itsDensityMatrix;
}

template class IrrepCD_Core<double>;
template class IrrepCD_Core<dcmplx>;
template class IrrepCD_HFPair<IrrepCD<double>>;
template class IrrepCD<double>;   // the FINITE leaf exists for double alone (no finite complex density)

// --- THE PERIODIC LEAF (both scalars; Step 3c-2b).  Note what is NOT here: no HF denials, no AO-face
// denial -- the leaf simply never inherits those capabilities, so nothing has to be denied (the R2.8
// smell the old IrrepCD<dcmplx> specializations carried is GONE).

// The complex overlap bypasses the (symmetric double) cache: MakeOverlap directly.  The REAL periodic
// block's S caches fine (theCache<double>), so it takes the core's cached accessor path.
template <class T> double PeriodicIrrepCD<T>::GetTotalCharge() const
{
    if constexpr (std::is_same_v<T,dcmplx>)
        return std::real(blazem::sum(this->itsDensityMatrix % blazem::trans(this->itsBasisSet->MakeOverlap())));
    else
        return IrrepCD_Core<T>::GetTotalCharge();
}

// NOT IMPLEMENTED (not "zero"): the periodic path is LDA-only; a silent rvec3_t(0,0,0) would hand a GGA
// a plausible wrong field (R1.4).
template <class T> rvec3_t PeriodicIrrepCD<T>::Gradient(const rvec3_t&) const
{
    throw std::logic_error("PeriodicIrrepCD::Gradient: grad(rho) is not wired for the periodic "
                           "density -- the periodic path is LDA-only.  Implement the gradient "
                           "contraction here, or ask the density through a gradient-capable face.");
}

// The run-typed energy contraction (Step 3c-2b): this REAL block's D against the complex run's dynamic
// term, through the term's Dynamic_CC_RealBlock face.  Same trace identity as the core's DM_Contract.
template <class T> double PeriodicIrrepCD<T>::DM_ContractE(const Dynamic_CC_RealBlock* v,
                                                           const tChargeDensity<dcmplx>* cd) const
{
    if constexpr (std::is_same_v<T,double>)
        return blazem::sum(this->itsDensityMatrix % blazem::trans(v->GetEMatrixR(this->itsBasisSet,this->itsSpin,cd)));
    else
        throw std::logic_error("PeriodicIrrepCD<dcmplx>::DM_ContractE: the run-typed contraction serves "
                               "REAL blocks inside a complex run; a complex block contracts natively");
}

template class IrrepCD_Fourier<PeriodicIrrepCD<double>>;
template class IrrepCD_Fourier<PeriodicIrrepCD<dcmplx>>;
template class PeriodicIrrepCD<double>;
template class PeriodicIrrepCD<dcmplx>;




rvec3_t GradientContraction(const vec_t<rvec3_t>& g, const vec_t<double>& v, const rsmat_t& m)
{
    // No UT coverage
    assert(v.size      ()==m.columns());
    assert(g.size      ()==m.columns());

    rvec3_t ret(0,0,0);
    for (unsigned int i=0; i<v.size(); i++)
        for (unsigned int j=0; j<v.size(); j++)
            ret+=m(i,j)*(g[i]*v[j]+v[i]*g[j]);
    return ret;
}

// rvec3_t GradientContraction(const Vector<Vec3>& g, const Vector<std::complex<double> >& v, const smat_t<std::complex<double> >& m)
// {
//     // No UT coverage
//     assert(v.size      ()==m.columns());
//     assert(g.size      ()==m.columns());

//     Vec3 ret(0,0,0);
//     for (unsigned int i=1; i<=v.size(); i++)
//         for (unsigned int j=1; j<=v.size(); j++)
//             ret+=m(i-1,j-1)*(g(i)*conj(v(j))+v(i)*conj(g(j)));
//     return real(ret);
// }

} //namespace
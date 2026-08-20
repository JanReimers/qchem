// File: BasisSet/Molecule/PG_Spherical/Imp/LatticeView.C  (see the interface unit for the design)
//
// Everything here is UNEXPORTED: consumers see only MakeSphericalLatticeView -> Real_BS.  The view
// holds its inner basis through the ABSTRACT faces alone (Molecule::Orbital_1E_IBS + LatticeSum1E),
// so it is engine-blind (MnD today, libCint tomorrow) and family-loose by construction.
module;
#include <cassert>
#include <cmath>       // std::sqrt (the column normalisation)
#include <complex>     // std::real/std::conj (the Hermitian packs)
#include <functional>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.BasisSet.Molecule.PG_Spherical.LatticeView;
import qchem.BasisSet;                        // Real_BS / tBasisSet<double>
import qchem.BasisSet.Orbital_1E_IBS;         // Real_OIBS faces (IDs, 1E integrals, VectorFunction)
import qchem.BasisSet.Molecule.IBS;           // Molecule::Orbital_1E_IBS (GetAoShells -- the inner shell source)
import qchem.BasisSet.Molecule.LatticeSum1E;  // the periodic capability the view answers
import qchem.BasisSet.Internal.BasisSetImp;   // the one-block Real_BS container
import qchem.Symmetry.Molecule.OperationRep;  // AoShell (+ ShellRep::Monomials soft capability)
import qchem.Math.Angular;                    // Math::SphericalShell / Monomial (the C_l source of truth)
import qchem.Structure;
import qchem.UnitCell;
import qchem.Blaze;
import qchem.Math;

namespace qchem::BasisSet::Molecule::PG_Spherical
{
using Symmetry::Molecule::AoShell;

//--------------------------------------------------------------------------------------------------
// The cart->sphere transform T (nCart x nSph), block-diagonal per shell.  A shell whose monomial
// degree is <=1 (s, p) passes through as identity -- a Cartesian s/p IS a spherical s/p (kept in the
// basis's own component order).  A degree-L>=2 Cartesian shell ((L+1)(L+2)/2 components) maps to its
// 2L+1 real solid harmonics: column coefficients c_t/N_c from Math::SphericalShell's RAW expansion
// (the code's single source of truth for the m-convention), then each column is normalised against
// the INNER shell's own overlap block -- so the normalisation convention is measured, not re-derived.
static rmat_t BuildCartToSphere(const std::vector<AoShell>& shells, const hmat_t<double>& S, size_t nCart)
{
    // Pass 1: the spherical column count.
    size_t nSph=0;
    std::vector<int> shellL(shells.size());
    for (size_t s=0;s<shells.size();s++)
    {
        const auto monos=shells[s].rep->Monomials();
        if (monos.empty())
            throw std::runtime_error("SphericalLatticeView: inner shell "+std::to_string(s)
                +" does not expose Cartesian monomials (ShellRep::Monomials) -- only a Cartesian-monomial "
                 "basis can be viewed spherically");
        const int L=monos.front().n+monos.front().l+monos.front().m;
        for (const auto& m : monos)
            if (m.n+m.l+m.m != L)
                throw std::runtime_error("SphericalLatticeView: mixed-degree shell (not a pure-L Cartesian shell)");
        shellL[s]=L;
        if (L<=1) nSph+=monos.size();
        else
        {
            if (size_t((L+1)*(L+2)/2)!=monos.size())
                throw std::runtime_error("SphericalLatticeView: degree-"+std::to_string(L)
                    +" shell with "+std::to_string(monos.size())+" components is not a full Cartesian shell");
            nSph+=size_t(2*L+1);
        }
    }

    // Pass 2: fill the blocks.
    rmat_t T(nCart,nSph,0.0);
    size_t col=0;
    for (size_t s=0;s<shells.size();s++)
    {
        const AoShell& sh=shells[s];
        const auto monos=sh.rep->Monomials();
        const size_t nc=monos.size(), off=sh.offset;
        const int L=shellL[s];
        if (L<=1)
        {
            for (size_t c=0;c<nc;c++) T(off+c,col+c)=1.0;
            col+=nc;
            continue;
        }
        const Math::HarmonicC2S c2s=Math::SphericalShell(L);      // 2L+1 harmonics, RAW relative weights
        for (const auto& harmonic : c2s)
        {
            rvec_t v(nc,0.0);
            for (const auto& t : harmonic)
            {
                size_t c=nc;                                       // find the term's component in THIS shell's order
                for (size_t i=0;i<nc;i++) if (monos[i]==t.p) { c=i; break; }
                if (c==nc) throw std::runtime_error("SphericalLatticeView: harmonic monomial missing from shell");
                v[c]=t.c/sh.norm[c];                               // AO_c = N_c (mono_c radial)  =>  mono_c radial = AO_c/N_c
            }
            double nrm2=0.0;                                       // S-block normalisation: nrm^2 = v^T S_block v
            for (size_t a=0;a<nc;a++)
                for (size_t b=0;b<nc;b++) nrm2+=v[a]*S(off+a,off+b)*v[b];
            assert(nrm2>0.0);
            const double inv=1.0/std::sqrt(nrm2);
            for (size_t c=0;c<nc;c++) T(off+c,col)=v[c]*inv;
            ++col;
        }
    }
    assert(col==nSph);
    return T;
}

//--------------------------------------------------------------------------------------------------
// The view block: Real_OIBS + GetAoShells + LatticeSum1E in the spherical span, everything by
// congruence/forward over the inner faces.  Hermitian/symmetric packing is explicit (upper triangle,
// real diagonal), mirroring the GPW evaluator's own convention.
class SphericalView_IBS
    : public virtual Molecule::Orbital_1E_IBS
    , public virtual Molecule::LatticeSum1E
{
public:
    SphericalView_IBS(std::shared_ptr<const Real_BS> holder, const Molecule::Orbital_1E_IBS* obs,
                      const Molecule::LatticeSum1E* lat, rmat_t T)
        : itsHolder(std::move(holder)), itsObs(obs), itsLat(lat), itsT(std::move(T))
        , itsTc(itsT.rows(),itsT.columns())
    {
        for (size_t i=0;i<itsT.rows();i++)
            for (size_t j=0;j<itsT.columns();j++) itsTc(i,j)=dcmplx(itsT(i,j),0.0);
        // T's NONZEROS, per spherical function (2026-08-20).  T is block-diagonal per shell and IDENTITY
        // for every s/p shell, so a dense trans(T)*v does ~nCart*nSph flops to compute what is, per output,
        // a 1-to-6 term sum.  On the GPW mesh path that dense product ran once per LATTICE IMAGE per mesh
        // point -- ~150x more often than the physics needs, and it was 10.6% of the whole run's cycles.
        // Ascending c, so the accumulation order matches the dense product's exactly (the skipped terms
        // contribute a hard 0.0): same values, bit for bit.
        itsTnz.resize(itsT.columns());
        for (size_t s=0;s<itsT.columns();s++)
            for (size_t c=0;c<itsT.rows();c++)
                if (itsT(c,s)!=0.0) itsTnz[s].push_back({c,itsT(c,s)});
    }

    // ---- IDs / symmetry / bookkeeping (forward; the ID suffix keeps every cache identity distinct) ----
    virtual std::string Name      () const override {return itsObs->Name()+"_sph";}
    virtual std::string BasisSetID() const override {return itsObs->BasisSetID()+"|sph";}
    virtual const Symmetry::Symmetry& GetSymmetry() const override {return itsObs->GetSymmetry();}
    virtual const sym_t& GetSymt  () const override {return itsObs->GetSymt();}
    virtual Irrep GetIrrep(const Spin& s) const override {return itsObs->GetIrrep(s);}
    virtual size_t GetNumFunctions() const override {return itsT.columns();}
    virtual std::ostream& Write(std::ostream& os) const override
    {
        os << "Spherical lattice view ("<<GetNumFunctions()<<" of "<<itsT.rows()<<" Cartesian) over: ";
        return itsObs->Write(os);
    }
    //! SALC over the view is doc/SphericalLatticePlan.md I4 (needs spherical AoShells); loud until then.
    virtual std::vector<AoShell> GetAoShells() const override
    { throw std::runtime_error("SphericalLatticeView: GetAoShells (SALC adaptation) not yet supported (plan I4)"); }

    // ---- point values (the GPW mesh path): v_sph = T^T v_cart ----
    virtual rvec_t operator()(const rvec3_t& r) const override
    {
        const rvec_t v=(*itsObs)(r);
        rvec_t out(itsT.columns(), 0.0);          // v_sph = T^T v_cart, over T's nonzeros only (see itsTnz)
        for (size_t s=0;s<out.size();s++)
        {
            double a=0.0;
            for (const auto& [c,t] : itsTnz[s]) a+=t*v[c];
            out[s]=a;
        }
        return out;
    }
    virtual rvec3vec_t Gradient(const rvec3_t& r) const override
    {
        const rvec3vec_t g=itsObs->Gradient(r);
        const size_t nC=itsT.rows(), nS=itsT.columns();
        rvec3vec_t out(nS, rvec3_t(0,0,0));
        for (size_t c=0;c<nC;c++)
            for (size_t s=0;s<nS;s++)
                if (itsT(c,s)!=0.0) out[s]=out[s]+itsT(c,s)*g[c];
        return out;
    }

    // ---- finite 1E (the home-only mode + the vet): congruence over the inner CACHED accessors ----
    virtual hmat_t<double> MakeOverlap() const override {return CongruenceR(itsObs->Overlap());}
    virtual hmat_t<double> MakeKinetic() const override {return CongruenceR(itsObs->Kinetic());}
    virtual hmat_t<double> MakeNuclear(const Structure* cl) const override {return CongruenceR(itsObs->Nuclear(cl));}

    // ---- LatticeSum1E: matrices/vectors by congruence, densities expanded, policy faces forwarded ----
    virtual chmat_t MakeOverlap(const cellphase_t& phase, const UnitCell& A) const override
    { return CongruenceC(itsLat->MakeOverlap(phase,A)); }
    virtual cvec_t  MakeOverlap(const cellphase_t& phase, const UnitCell& A, const GaussianFunction& g) const override
    { cvec_t b=itsLat->MakeOverlap(phase,A,g); cvec_t out=blazem::trans(itsTc)*b; return out; }
    virtual cvec_t  MakeOverlap(const GaussianFunction& g) const override
    { cvec_t b=itsLat->MakeOverlap(g); cvec_t out=blazem::trans(itsTc)*b; return out; }
    // ★ THE POINT OF THE POINT-SET SEAM.  The transform is LINEAR, so
    //       Sum_R phase_R T^T v_R  ==  T^T Sum_R phase_R v_R,
    // and the right-hand side applies T ONCE PER POINT where operator() applied it once per IMAGE (~150x
    // on the MnO benchmark cell).  This override is only expressible because the inner basis does the
    // image sum: with the sum in the periodic CALLER, the caller holds an abstract basis and cannot know
    // a transform is hiding inside it (doc/OpenWork.md Step 3).  Same nonzero-only contraction as
    // operator() -- T is block-diagonal per shell and IDENTITY for every s/p shell.
    virtual void BlochPointValues(const rvec3vec_t& pts, const cellphase_t& phase, const UnitCell& A,
                                  mat_t<dcmplx>& Phi) const override
    {
        mat_t<dcmplx> C;                                  // the inner CARTESIAN Bloch table
        itsLat->BlochPointValues(pts,phase,A,C);
        const size_t npts=C.rows(), nS=itsT.columns();
        Phi=mat_t<dcmplx>(npts,nS,dcmplx(0.0));
        for (size_t q=0;q<npts;q++)
            for (size_t s=0;s<nS;s++)
            {
                dcmplx a(0.0);
                for (const auto& [c,t] : itsTnz[s]) a+=t*C(q,c);
                Phi(q,s)=a;
            }
    }
    virtual chmat_t MakeKinetic(const cellphase_t& phase, const UnitCell& A) const override
    { return CongruenceC(itsLat->MakeKinetic(phase,A)); }
    virtual chmat_t MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const override
    { return CongruenceC(itsLat->MakeNuclear(phase,A,cl)); }
    virtual chmat_t MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                                      const std::function<GaussianFunction(int)>& opForZ) const override
    { return CongruenceC(itsLat->MakeLocalGaussian(phase,A,cl,opForZ)); }
    virtual chmat_t MakeLocalGaussian(const Structure* cl,
                                      const std::function<GaussianFunction(int)>& opForZ) const override
    { return CongruenceC(itsLat->MakeLocalGaussian(cl,opForZ)); }
    virtual double MaxExponent    () const override {return itsLat->MaxExponent();}
    virtual double MinExponent    () const override {return itsLat->MinExponent();}
    virtual double RelCutoffSafety() const override {return itsLat->RelCutoffSafety();}
    virtual std::vector<rvec_t> CollocateDensity(const chmat_t& D, const cellphase_t& phase, const UnitCell& A,
                                                 const std::vector<ivec3_t>& N_L,
                                                 const std::vector<double>& ecut_L,
                                                 double relFieldSharp=-1.0) const override
    { return itsLat->CollocateDensity(Expand(D),phase,A,N_L,ecut_L,relFieldSharp); }
    virtual chmat_t IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                                       const std::vector<ivec3_t>& N_L,
                                       const std::vector<double>& ecut_L, double absRelCutoff=0.0,
                                       const chmat_t* screenD=nullptr, double fieldSharpness=0.0,
                                       double relFieldSharp=-1.0,
                                       const std::vector<size_t>* pairLevels=nullptr) const override
    {
        chmat_t Dexp; const chmat_t* sd=nullptr;
        if (screenD) { Dexp=Expand(*screenD); sd=&Dexp; }
        // pairLevels is TOKEN-PASSED: produced by our StaticFieldPairLevels forward, so it is already
        // inner(Cartesian)-pair-indexed -- the two forwards stay a consistent pair by construction.
        return CongruenceC(itsLat->IntegratePotential(V_L,phase,A,N_L,ecut_L,absRelCutoff,sd,
                                                      fieldSharpness,relFieldSharp,pairLevels));
    }
    virtual std::vector<size_t> StaticFieldPairLevels(const std::vector<double>& ecut_L,
                                                      double beta, double lnEps) const override
    { return itsLat->StaticFieldPairLevels(ecut_L,beta,lnEps); }
    virtual void EmitLatticeSumReport(const UnitCell& A) const override {itsLat->EmitLatticeSumReport(A);}
    virtual void ReleaseStreams(const std::vector<ivec3_t>& N_L, const std::vector<double>& ecut_L) const override
    { itsLat->ReleaseStreams(N_L,ecut_L); }
    virtual size_t SetStreamSymmetryOps(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops,
                                        const UnitCell& A, const rvec3_t& kFrac=rvec3_t(0,0,0)) const override
    { return itsLat->SetStreamSymmetryOps(ops,A,kFrac); }   // the fold acts on the INNER Cartesian streams
    virtual size_t StreamFoldOrder() const override {return itsLat->StreamFoldOrder();}

private:
    // T^T M T, packed (upper triangle, real diagonal for the complex case).
    hmat_t<double> CongruenceR(const hmat_t<double>& M) const
    {
        const size_t nS=itsT.columns();
        rmat_t A=blazem::trans(itsT)*rmat_t(M)*itsT;
        hmat_t<double> H(nS);
        for (size_t i=0;i<nS;i++) for (size_t j=i;j<nS;j++) H(i,j)=0.5*(A(i,j)+A(j,i));
        return H;
    }
    chmat_t CongruenceC(const chmat_t& M) const
    {
        const size_t nS=itsTc.columns();
        mat_t<dcmplx> A=blazem::trans(itsTc)*mat_t<dcmplx>(M)*itsTc;
        chmat_t H=blazem::zeroH<dcmplx>(nS);
        for (size_t i=0;i<nS;i++)
        {
            H(i,i)=dcmplx(std::real(A(i,i)),0.0);
            for (size_t j=i+1;j<nS;j++) H(i,j)=0.5*(A(i,j)+std::conj(A(j,i)));
        }
        return H;
    }
    chmat_t Expand(const chmat_t& Dsph) const   // D_cart = T D_sph T^T
    {
        const size_t nC=itsTc.rows();
        mat_t<dcmplx> A=itsTc*mat_t<dcmplx>(Dsph)*blazem::trans(itsTc);
        chmat_t H=blazem::zeroH<dcmplx>(nC);
        for (size_t i=0;i<nC;i++)
        {
            H(i,i)=dcmplx(std::real(A(i,i)),0.0);
            for (size_t j=i+1;j<nC;j++) H(i,j)=0.5*(A(i,j)+std::conj(A(j,i)));
        }
        return H;
    }

    std::shared_ptr<const Real_BS> itsHolder;   //!< keeps the inner basis (and its evaluator) alive
    const Molecule::Orbital_1E_IBS* itsObs;     //!< the inner orbital faces (abstract)
    const Molecule::LatticeSum1E*   itsLat;     //!< the inner periodic capability (abstract)
    rmat_t        itsT;                         //!< cart->sphere, nCart x nSph
    mat_t<dcmplx> itsTc;                        //!< the same T, complex, for the chmat congruences
    //! T's nonzeros per spherical function: (cartesian index, coefficient), ascending index.  The
    //! POINTWISE path only -- the matrix congruences keep the dense form, where blaze's blocked kernels
    //! are the right tool and the call count is O(1) per run rather than O(points x images).
    std::vector<std::vector<std::pair<size_t,double>>> itsTnz;
};

//--------------------------------------------------------------------------------------------------
// The one-block Real_BS container (mirrors PG_Cart::BasisSet's shape).
class ViewBS
    : public virtual ::qchem::BasisSet::tBasisSet<double>
    , public ::qchem::BasisSet::BasisSetImp<double>
{
public:
    ViewBS() {}
    virtual void Insert(obs_t* bs) { ::qchem::BasisSet::BasisSetImp<double>::Insert(bs); }
};

std::shared_ptr<const Real_BS> MakeSphericalLatticeView(std::shared_ptr<const Real_BS> cart)
{
    const Molecule::Orbital_1E_IBS* obs=nullptr;
    for (auto ibs : const_cast<Real_BS&>(*cart).Iterate<Molecule::Orbital_1E_IBS>()) { obs=ibs; break; }
    if (!obs) throw std::runtime_error("MakeSphericalLatticeView: no Molecule::Orbital_1E_IBS block in the wrapped basis");
    const auto* lat=dynamic_cast<const Molecule::LatticeSum1E*>(obs);
    if (!lat) throw std::runtime_error("MakeSphericalLatticeView: the wrapped orbital block has no LatticeSum1E capability");
    const rmat_t T=BuildCartToSphere(obs->GetAoShells(), obs->Overlap(), obs->GetNumFunctions());
    auto* bs=new ViewBS;
    bs->Insert(new SphericalView_IBS(cart, obs, lat, T));
    return std::shared_ptr<const Real_BS>(bs);
}

} //namespace

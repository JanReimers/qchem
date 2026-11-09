// File: BasisSet/Gaussian/Evaluators/PG_LibCint/Imp/Evaluator.C
//
// libcint packing + matrix assembly for the Cartesian PG basis.  See Evaluator.C for the why; this file
// is the how:
//   1. Regroup PGData's flattened (radial x polarization) components back into SHELLS (a run of components
//      sharing one radial and one total-L).  Each shell becomes one libcint contracted Cartesian basis
//      function (NCTR=1) with the radial's exponents + (normalization-folded) contraction coeffs.
//   2. Per shell, build the libcint-cart -> PG-order component permutation (lc2pg) by matching (lx,ly,lz).
//   3. Renormalize every component by 1/sqrt(native self-overlap) so the delivered matrices use PG's
//      unit-self-overlap convention (== PGData::ns, but computed from libcint to stay convention-free).
// Matrices are then assembled by looping libcint shell tuples (the library's natural grain), scattering
// each dense block into the PG-ordered result with the permutation + scales.
module;
#include <vector>
#include <cmath>
#include <cassert>
#include <ostream>
#include <array>
#include <stdexcept>
#include <string>
extern "C" {
#include "cint_funcs.h"   // pulls cint.h: slot constants, FINT/CACHE_SIZE_T, the int*_cart functions
// cint_funcs.h omits prototypes for the 3-centre integrals (the symbols are in libcint); declare them.
CINTIntegralFunction int3c1e_cart;   // 3-centre overlap  <ab|c>
CINTIntegralFunction int3c2e_cart;   // 3-centre Coulomb  (ab|c)
CINTIntegralFunction int2c2e_cart;   // 2-centre Coulomb  (a|b)
CINTIntegralFunction int3c1e_ip1_cart;   // <∇a b|c>: the local PP's short-part gradient
double CINTcommon_fac_sp(FINT l);      // the s/p factor libcint folds into every integral (g1e.h, not in the public header)
}
module qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;
import qchem.Structure;
import qchem.Types;
import qchem.Math.SolidHarmonics;   // Math::CartGaussian
import qchem.Math.Angular;          // Math::CartTerm / Monomial (the charge Gaussian's one term)

namespace qchem::BasisSet::Gaussian::Evaluators::PG_LibCint
{
using PG_Cart_MnD::PGData;
using PG_Cart_MnD::Polarization;
using PG_Cart_MnD::GaussianRF;
using IntFn = CINTIntegralFunction*;   // (out,dims,shls,atm,natm,bas,nbas,env,opt,cache)

static int ncart(int l) {return (l+1)*(l+2)/2;}

// --- internal layout ------------------------------------------------------------------------------
struct AtomRec { int Z; double x,y,z; };
struct ShellDesc
{
    int                 atom;   // index into the (shared) atom table
    int                 l;
    std::vector<double> exp;    // primitive exponents
    std::vector<double> coef;   // PG normalization-folded contraction coefficients
    std::vector<size_t> lc2pg;  // libcint component -> global output index (Cartesian: PG order; spherical:
                                // sequential = libcint's native real-spherical order)
    std::vector<double> c2s;    // spherical only: (2l+1) x ncart row-major, libcint's cart->sph map (grid-eval)
};

struct NR_Evaluator::Imp
{
    bool                   spherical=false;
    std::vector<AtomRec>   atoms;
    std::vector<ShellDesc> shells;
    rvec_t                 scale;   // per global component: 1/sqrt(native self-overlap)
    size_t                 N=0;     // number of delivered components ((l+1)(l+2)/2 cart, 2l+1 sph, per shell)

    int dim(int l) const {return spherical ? 2*l+1 : (l+1)*(l+2)/2;}   // components per shell

    // Build atm/bas/env for the shared atom table + a concatenation of shell lists (slots).  Shell ids in
    // the returned bas run in the order the ShellDesc pointers are given.
    static void BuildMol(const std::vector<AtomRec>& atoms, const std::vector<const ShellDesc*>& sh,
                         std::vector<int>& atm, std::vector<int>& bas, std::vector<double>& env)
    {
        atm.assign(atoms.size()*ATM_SLOTS, 0);
        env.assign(PTR_ENV_START, 0.0);
        for (size_t i=0;i<atoms.size();++i)
        {
            atm[CHARGE_OF+ATM_SLOTS*i]=atoms[i].Z;
            atm[PTR_COORD+ATM_SLOTS*i]=(int)env.size();
            env.push_back(atoms[i].x); env.push_back(atoms[i].y); env.push_back(atoms[i].z);
        }
        bas.assign(sh.size()*BAS_SLOTS, 0);
        for (size_t s=0;s<sh.size();++s)
        {
            bas[ATOM_OF +BAS_SLOTS*s]=sh[s]->atom;
            bas[ANG_OF  +BAS_SLOTS*s]=sh[s]->l;
            bas[NPRIM_OF+BAS_SLOTS*s]=(int)sh[s]->exp.size();
            bas[NCTR_OF +BAS_SLOTS*s]=1;
            bas[PTR_EXP +BAS_SLOTS*s]=(int)env.size();
            for (double e:sh[s]->exp ) env.push_back(e);
            bas[PTR_COEFF+BAS_SLOTS*s]=(int)env.size();
            for (double c:sh[s]->coef) env.push_back(c);
        }
    }
};

// libcint Cartesian component order for momentum l (CINTcart_comp: lx desc, ly desc), as Polarizations.
static std::vector<Polarization> LibcintCartOrder(int l)
{
    std::vector<Polarization> v;
    for (int lx=l; lx>=0; --lx)
        for (int ly=l-lx; ly>=0; --ly)
            v.push_back(Polarization(lx, ly, l-lx-ly));   // (n,l,m) = (x,y,z) powers
    return v;
}

NR_Evaluator::NR_Evaluator() = default;

NR_Evaluator::NR_Evaluator(const PGData& data, const Structure* cl)
{
    static_cast<PGData&>(*this) = data;             // own the component layout + order (size/Norm/ns)
    Init(cl);
}

// Build the libcint shell + atom tables from this evaluator's (already-populated) PGData and the structure.
void NR_Evaluator::Init(const Structure* cl, bool spherical)
{
    itsImp.reset(new Imp);
    Imp& m = *itsImp;
    m.spherical = spherical;
    size_t Ncart = PGData::size();                  // # Cartesian components flattened in the PGData

    // Atom table (shell centres + nuclear charges) from the structure.
    for (const Site& a : cl->Sites()) m.atoms.push_back(AtomRec{a.Z, a.R.x, a.R.y, a.R.z});
    auto findAtom = [&](const rvec3_t& R)->int
    {
        for (size_t i=0;i<m.atoms.size();++i)
            if (std::fabs(m.atoms[i].x-R.x)+std::fabs(m.atoms[i].y-R.y)+std::fabs(m.atoms[i].z-R.z) < 1e-8)
                return (int)i;
        assert(false && "shell centre does not match a nucleus"); return 0;
    };

    // Regroup the flattened Cartesian components into shells: a new shell begins when the radial pointer OR
    // the total L changes (PGData lays out block-by-block, MakePolarizations groups by L, so runs are
    // contiguous).  Packing into libcint is the SAME (Cartesian shell defs) for both modes; only the
    // delivered component count + order differ (Cartesian -> PG permutation; spherical -> 2l+1 sequential).
    size_t gout=0;
    for (size_t i=0;i<Ncart;)
    {
        const GaussianRF* r = radials[i];
        int L = pols[i].GetTotalL();
        size_t j=i;
        std::vector<size_t> gidx;                   // global Cartesian component indices (PG order)
        while (j<Ncart && radials[j]==r && pols[j].GetTotalL()==L) gidx.push_back(j++);
        assert((int)gidx.size()==ncart(L) && "shell is not a complete Cartesian set");

        ShellDesc sh;
        sh.atom = findAtom(r->GetCenter());
        sh.l    = L;
        { rvec_t es=r->GetExponents(); const rvec_t& cs=r->GetCoeffs();
          for (size_t k=0;k<es.size();++k){ sh.exp.push_back(es[k]); sh.coef.push_back(cs[k]); } }
        if (spherical)
        {
            for (int p=0;p<m.dim(L);++p) sh.lc2pg.push_back(gout++);    // libcint-native sph order, sequential
            // The cart->sph map is libcint's own: CINTc2s_bra_sph of the identity (one Cartesian unit vector per
            // "ket" row).  s and p return the INPUT pointer (no copy), so read through the return value.
            const int nc=ncart(L), ns=2*L+1;
            std::vector<double> id((size_t)nc*nc, 0.0), out((size_t)nc*ns, 0.0);
            for (int c=0;c<nc;++c) id[(size_t)c*nc+c]=1.0;
            const double* t=CINTc2s_bra_sph(out.data(), nc, id.data(), L);
            sh.c2s.assign((size_t)ns*nc, 0.0);
            for (int c=0;c<nc;++c) for (int p=0;p<ns;++p) sh.c2s[(size_t)p*nc+c]=t[(size_t)c*ns+p];
        }
        else
            for (const Polarization& lp:LibcintCartOrder(L))           // libcint-cart order -> PG global index
            {
                size_t hit=gidx.size();
                for (size_t k=0;k<gidx.size();++k) if (pols[gidx[k]]==lp){ hit=k; break; }
                assert(hit<gidx.size() && "missing Cartesian component"); sh.lc2pg.push_back(gidx[hit]);
            }
        m.shells.push_back(std::move(sh));
        i=j;
    }
    m.N = spherical ? gout : Ncart;

    // Per-component scale = 1/sqrt(native self-overlap) (libcint does not unit-normalize Cartesians; the
    // _sph harmonics are already normalized, so this is ~1 there, but compute it anyway -- convention-free).
    m.scale.resize(m.N);
    std::vector<const ShellDesc*> self; for (auto& s:m.shells) self.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(m.atoms, self, atm,bas,env);
    int natm=(int)m.atoms.size(), nbas=(int)m.shells.size();
    IntFn ovlp = spherical ? int1e_ovlp_sph : int1e_ovlp_cart;
    std::vector<double> buf;
    for (size_t s=0;s<m.shells.size();++s)
    {
        int d=m.dim(m.shells[s].l); buf.assign((size_t)d*d, 0.0);
        int shls[2]={(int)s,(int)s};
        ovlp(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
        for (int p=0;p<d;++p) m.scale[m.shells[s].lc2pg[p]] = 1.0/std::sqrt(buf[p+d*p]);
    }
}

NR_Evaluator::~NR_Evaluator() = default;

// --- spherical grid-eval: phi_p(r) = scale_p * fac_sp(l) * sum_c c2s[p][c] x^lx y^ly z^lz * sum_k coef_k exp(-a_k r^2) ----
// The factors are exactly libcint's (its integrals carry fac_sp and c2s), so the delivered scale normalizes these to one.
static double IPow(double x, int n) {double r=1.0; for (int i=0;i<n;++i) r*=x; return r;}

rvec_t NR_Evaluator::SphericalValues(const rvec3_t& r) const
{
    const Imp& m=*itsImp;
    assert(m.spherical && "SphericalValues on a Cartesian libcint evaluator");
    rvec_t ret(m.N);
    std::vector<double> mono;
    for (const ShellDesc& sh:m.shells)
    {
        const AtomRec& A=m.atoms[sh.atom];
        const double dx=r.x-A.x, dy=r.y-A.y, dz=r.z-A.z, r2=dx*dx+dy*dy+dz*dz;
        double R=0.0;
        for (size_t k=0;k<sh.exp.size();++k) R+=sh.coef[k]*std::exp(-sh.exp[k]*r2);
        R*=CINTcommon_fac_sp(sh.l);
        const int l=sh.l, nc=ncart(l), ns=2*l+1;
        mono.assign(nc, 0.0);
        int c=0;
        for (int lx=l; lx>=0; --lx)                                   // libcint's Cartesian order (lx desc, ly desc)
            for (int ly=l-lx; ly>=0; --ly) mono[c++]=IPow(dx,lx)*IPow(dy,ly)*IPow(dz,l-lx-ly);
        for (int p=0;p<ns;++p)
        {
            double P=0.0;
            for (int q=0;q<nc;++q) P+=sh.c2s[(size_t)p*nc+q]*mono[q];
            const size_t g=sh.lc2pg[p];
            ret[g]=m.scale[g]*P*R;
        }
    }
    return ret;
}

rvec3vec_t NR_Evaluator::SphericalGradient(const rvec3_t& r) const
{
    const Imp& m=*itsImp;
    assert(m.spherical && "SphericalGradient on a Cartesian libcint evaluator");
    rvec3vec_t ret(m.N);
    std::vector<double> mono, gx, gy, gz;
    for (const ShellDesc& sh:m.shells)
    {
        const AtomRec& A=m.atoms[sh.atom];
        const double dx=r.x-A.x, dy=r.y-A.y, dz=r.z-A.z, r2=dx*dx+dy*dy+dz*dz;
        double R=0.0, dR=0.0;                                         // grad R = dR * (dx,dy,dz)
        for (size_t k=0;k<sh.exp.size();++k)
        {
            const double e=sh.coef[k]*std::exp(-sh.exp[k]*r2);
            R+=e; dR+=-2.0*sh.exp[k]*e;
        }
        const double f=CINTcommon_fac_sp(sh.l); R*=f; dR*=f;
        const int l=sh.l, nc=ncart(l), ns=2*l+1;
        mono.assign(nc,0.0); gx.assign(nc,0.0); gy.assign(nc,0.0); gz.assign(nc,0.0);
        int c=0;
        for (int lx=l; lx>=0; --lx)
            for (int ly=l-lx; ly>=0; --ly, ++c)
            {
                const int lz=l-lx-ly;
                const double px=IPow(dx,lx), py=IPow(dy,ly), pz=IPow(dz,lz);
                mono[c]=px*py*pz;
                gx[c]=(lx>0 ? lx*IPow(dx,lx-1) : 0.0)*py*pz;
                gy[c]=px*(ly>0 ? ly*IPow(dy,ly-1) : 0.0)*pz;
                gz[c]=px*py*(lz>0 ? lz*IPow(dz,lz-1) : 0.0);
            }
        for (int p=0;p<ns;++p)
        {
            double P=0.0, Px=0.0, Py=0.0, Pz=0.0;
            for (int q=0;q<nc;++q)
            {
                const double w=sh.c2s[(size_t)p*nc+q];
                P+=w*mono[q]; Px+=w*gx[q]; Py+=w*gy[q]; Pz+=w*gz[q];
            }
            const size_t g=sh.lc2pg[p];
            const double s=m.scale[g];
            ret[g]=rvec3_t(s*(Px*R+P*dR*dx), s*(Py*R+P*dR*dy), s*(Pz*R+P*dR*dz));
        }
    }
    return ret;
}

size_t NR_Evaluator::size() const {return itsImp ? itsImp->N     : PGData::size();}
rvec_t NR_Evaluator::Norm() const {return itsImp ? itsImp->scale : ns;}

std::ostream& NR_Evaluator::Write(std::ostream& os) const
{
    return os << "PG_LibCint::NR_Evaluator[" << size() << "]";
}

// --- 1E: loop shell pairs, scatter each block into the symmetric PG-ordered matrix ----------------
rsmat_t NR_Evaluator::Build1(int which) const
{
    const Imp& m=*itsImp;
    rsmat_t M(m.N);
    std::vector<const ShellDesc*> self; for (auto& s:m.shells) self.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(m.atoms, self, atm,bas,env);
    int natm=(int)m.atoms.size(), nbas=(int)m.shells.size();
    IntFn fn   = m.spherical
        ? (which==0 ? int1e_ovlp_sph  : which==1 ? int1e_kin_sph  : int1e_nuc_sph)
        : (which==0 ? int1e_ovlp_cart : which==1 ? int1e_kin_cart : int1e_nuc_cart);
    double pre = which==1 ? 2.0 : 1.0;   // KineticMatrix is <-nabla^2> = 2 x libcint's T = -1/2 nabla^2
    std::vector<double> buf;
    for (size_t si=0;si<m.shells.size();++si)
        for (size_t sj=si;sj<m.shells.size();++sj)
        {
            int di=m.dim(m.shells[si].l), dj=m.dim(m.shells[sj].l);
            buf.assign((size_t)di*dj, 0.0);
            int shls[2]={(int)si,(int)sj};
            fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
            for (int p=0;p<di;++p){ size_t gi=m.shells[si].lc2pg[p];
                for (int q=0;q<dj;++q){ size_t gj=m.shells[sj].lc2pg[q];
                    M(gi,gj) = pre * buf[p+di*q] * m.scale[gi]*m.scale[gj]; } }
        }
    return M;
}
// --- geometry gradients of the 1E matrices --------------------------------------------------------
// libcint's ip integrals differentiate the BRA function with respect to the electron coordinate: X(i,j) = <d_c i|O|j>.  A function centred
// on A has d/dA = -d/dr, so for a centre-only operator dM_ij/dR_A = -(X(i,j)[i on A] + X(j,i)[j on A]).  The nuclear operator also moves:
// with Y(i,j) = <d_c i| 1/|r - R_C| |j> (int1e_iprinv, its origin on C), translation invariance of the three-body integral gives
// d<i|1/r_C|j>/dR_C = Y(i,j) + Y(j,i), so dV/dR_C picks up -Z_C (Y + Y^T).
NR_Evaluator::GradientMatrices NR_Evaluator::BuildGradient(int which) const
{
    const Imp& m=*itsImp;
    if (m.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the 1E geometry gradients are Cartesian-only (no client asks for spherical yet)");
    std::vector<const ShellDesc*> self; for (auto& s:m.shells) self.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(m.atoms, self, atm,bas,env);
    const int natm=(int)m.atoms.size(), nbas=(int)m.shells.size();
    const double pre = which==1 ? 2.0 : 1.0;   // <-nabla^2> = 2 x libcint's T
    std::vector<int> compAtom(m.N, -1);
    for (const auto& sh : m.shells) for (size_t g : sh.lc2pg) compAtom[g]=sh.atom;
    // X_c(i,j) = <d_c i| O |j> over every ordered shell pair (not symmetric), scaled into the delivered normalization
    auto ip=[&](IntFn fn)
    {
        std::array<rmat_t,3> X; for (auto& x : X) x=rmat_t(m.N, m.N, 0.0);
        std::vector<double> buf;
        for (size_t si=0;si<m.shells.size();++si)
            for (size_t sj=0;sj<m.shells.size();++sj)
            {
                const int di=m.dim(m.shells[si].l), dj=m.dim(m.shells[sj].l);
                buf.assign((size_t)di*dj*3, 0.0);
                int shls[2]={(int)si,(int)sj};
                fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
                for (int c=0;c<3;++c)
                    for (int p=0;p<di;++p){ const size_t gi=m.shells[si].lc2pg[p];
                        for (int q=0;q<dj;++q){ const size_t gj=m.shells[sj].lc2pg[q];
                            X[c](gi,gj) = pre * buf[p+di*(q+dj*c)] * m.scale[gi]*m.scale[gj]; } }
            }
        return X;
    };
    GradientMatrices G(m.atoms.size());
    for (auto& a : G) for (auto& x : a) x=rmat_t(m.N, m.N, 0.0);
    IntFn fn = which==0 ? int1e_ipovlp_cart : which==1 ? int1e_ipkin_cart : int1e_ipnuc_cart;
    const std::array<rmat_t,3> X=ip(fn);
    for (int c=0;c<3;++c)
        for (size_t i=0;i<m.N;++i) for (size_t j=0;j<m.N;++j)
        {
            G[compAtom[i]][c](i,j)-=X[c](i,j);
            G[compAtom[j]][c](i,j)-=X[c](j,i);
        }
    if (which==2)   // the nuclei move too: one int1e_iprinv per nucleus, its origin on it
        for (size_t C=0;C<m.atoms.size();++C)
        {
            env[PTR_RINV_ORIG]=m.atoms[C].x; env[PTR_RINV_ORIG+1]=m.atoms[C].y; env[PTR_RINV_ORIG+2]=m.atoms[C].z;
            const std::array<rmat_t,3> Y=ip(int1e_iprinv_cart);
            const double Z=m.atoms[C].Z;
            for (int c=0;c<3;++c)
                for (size_t i=0;i<m.N;++i) for (size_t j=0;j<m.N;++j) G[C][c](i,j)-=Z*(Y[c](i,j)+Y[c](j,i));
        }
    return G;
}
NR_Evaluator::GradientMatrices NR_Evaluator::OverlapGradient() const {return BuildGradient(0);}
NR_Evaluator::GradientMatrices NR_Evaluator::KineticGradient() const {return BuildGradient(1);}
NR_Evaluator::GradientMatrices NR_Evaluator::NuclearGradient() const {return BuildGradient(2);}

rsmat_t NR_Evaluator::OverlapMatrix()             const {return Build1(0);}
rsmat_t NR_Evaluator::KineticMatrix()             const {return Build1(1);}
rsmat_t NR_Evaluator::NuclearMatrix(const Structure*) const {return Build1(2);} // atm already carries the nuclei

// --- 3-centre <ab|c>: one symmetric (ia,ib) block per fit component ic ----------------------------
Projector3<double> NR_Evaluator::Build3(const NR_Evaluator& fit, int which) const
{
    const Imp& a=*itsImp; const Imp& c=*fit.itsImp;
    assert(!a.spherical && "3-centre (DFT) is Cartesian-only for the libcint evaluator");
    std::vector<const ShellDesc*> sh; for (auto& s:a.shells) sh.push_back(&s);
    size_t cOff=sh.size();           for (auto& s:c.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(a.atoms, sh, atm,bas,env);
    int natm=(int)a.atoms.size(), nbas=(int)sh.size();
    IntFn fn = which==0 ? int3c1e_cart : int3c2e_cart;   // overlap : Coulomb

    Projector3<double> s3; s3.dense.reserve(c.N); for (size_t i=0;i<c.N;++i) s3.dense.push_back(rsmat_t(a.N));
    std::vector<double> buf;
    for (size_t sk=0;sk<c.shells.size();++sk)
        for (size_t si=0;si<a.shells.size();++si)
            for (size_t sj=si;sj<a.shells.size();++sj)
            {
                int di=a.dim(a.shells[si].l), dj=a.dim(a.shells[sj].l), dk=c.dim(c.shells[sk].l);
                buf.assign((size_t)di*dj*dk, 0.0);
                int shls[3]={(int)si,(int)sj,(int)(cOff+sk)};
                fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
                for (int p=0;p<di;++p){ size_t gia=a.shells[si].lc2pg[p];
                    for (int q=0;q<dj;++q){ size_t gib=a.shells[sj].lc2pg[q];
                        for (int rr=0;rr<dk;++rr){ size_t gic=c.shells[sk].lc2pg[rr];
                            s3.dense[gic](gia,gib) = buf[p+di*(q+dj*rr)] * a.scale[gia]*a.scale[gib]*c.scale[gic]; } } }
            }
    return s3;
}
// --- the Dunlap force's integrals, contracted -------------------------------------------------------
// ip1 differentiates the FIRST orbital function and ip2 the fit function, with respect to the electron coordinate; a centre derivative is
// minus that.  Over ORDERED orbital pairs (i, j), the ip1 sum is the first function's derivative; the second function's is the same sum
// relabelled (D and (ij|k) are symmetric in i, j), hence the factor 2.  ip2 runs over the same ordered pairs once.
std::vector<rvec3_t> NR_Evaluator::RepulsionThreeC_Gradient(const NR_Evaluator& fit, const rmat_t& D, const rvec_t& cf) const
{
    const Imp& a=*itsImp; const Imp& c=*fit.itsImp;
    if (a.spherical || c.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the 3-centre gradient is Cartesian-only");
    if (D.rows()!=a.N || D.columns()!=a.N || cf.size()!=c.N)
        throw std::invalid_argument("PG_LibCint::RepulsionThreeC_Gradient: D is "+std::to_string(D.rows())+"x"+std::to_string(D.columns())+" and c has "
                                    +std::to_string(cf.size())+" for "+std::to_string(a.N)+" orbital and "+std::to_string(c.N)+" fit functions");
    std::vector<const ShellDesc*> sh; for (auto& s:a.shells) sh.push_back(&s);
    const size_t cOff=sh.size();      for (auto& s:c.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(a.atoms, sh, atm,bas,env);
    const int natm=(int)a.atoms.size(), nbas=(int)sh.size();
    std::vector<rvec3_t> g(a.atoms.size(), rvec3_t(0,0,0));
    std::vector<double> buf;
    for (size_t sk=0;sk<c.shells.size();++sk)
        for (size_t si=0;si<a.shells.size();++si)
            for (size_t sj=0;sj<a.shells.size();++sj)
            {
                const int di=a.dim(a.shells[si].l), dj=a.dim(a.shells[sj].l), dk=c.dim(c.shells[sk].l);
                int shls[3]={(int)si,(int)sj,(int)(cOff+sk)};
                for (int which=0; which<2; ++which)
                {
                    buf.assign((size_t)di*dj*dk*3, 0.0);
                    (which==0 ? int3c2e_ip1_cart : int3c2e_ip2_cart)(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
                    const int on = which==0 ? a.shells[si].atom : c.shells[sk].atom;
                    double acc[3]={0,0,0};
                    for (int x=0;x<3;++x)
                        for (int p=0;p<di;++p){ const size_t gi=a.shells[si].lc2pg[p];
                            for (int q=0;q<dj;++q){ const size_t gj=a.shells[sj].lc2pg[q]; const double dij=D(gi,gj)*a.scale[gi]*a.scale[gj];
                                for (int r=0;r<dk;++r){ const size_t gk=c.shells[sk].lc2pg[r];
                                    acc[x]+=dij*cf[gk]*c.scale[gk]*buf[p+di*(q+dj*(r+dk*x))]; } } }
                    g[on]-=(which==0 ? 2.0 : 1.0)*rvec3_t(acc[0], acc[1], acc[2]);
                }
            }
    return g;
}

std::vector<rvec3_t> NR_Evaluator::RepulsionTwoC_Gradient(const rvec_t& cf) const
{
    const Imp& c=*itsImp;
    if (cf.size()!=c.N) throw std::invalid_argument("PG_LibCint::RepulsionTwoC_Gradient: c has "+std::to_string(cf.size())+" for "+std::to_string(c.N)+" fit functions");
    rsmat_t A(c.N);
    for (size_t k=0;k<c.N;++k) for (size_t l=k;l<c.N;++l) A(k,l)=cf[k]*cf[l];
    return RepulsionTwoC_Gradient(A);
}

std::vector<rvec3_t> NR_Evaluator::RepulsionTwoC_Gradient(const rsmat_t& A) const
{
    const Imp& c=*itsImp;
    if (c.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the 2-centre gradient is Cartesian-only");
    if (A.rows()!=c.N) throw std::invalid_argument("PG_LibCint::RepulsionTwoC_Gradient: the adjoint is "+std::to_string(A.rows())+" square for "+std::to_string(c.N)+" fit functions");
    std::vector<const ShellDesc*> sh; for (auto& s:c.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(c.atoms, sh, atm,bas,env);
    const int natm=(int)c.atoms.size(), nbas=(int)sh.size();
    std::vector<rvec3_t> g(c.atoms.size(), rvec3_t(0,0,0));
    std::vector<double> buf;
    for (size_t sk=0;sk<c.shells.size();++sk)
        for (size_t sl=0;sl<c.shells.size();++sl)
        {
            const int dk=c.dim(c.shells[sk].l), dl=c.dim(c.shells[sl].l);
            buf.assign((size_t)dk*dl*3, 0.0);
            int shls[2]={(int)sk,(int)sl};
            int2c2e_ip1_cart(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
            double acc[3]={0,0,0};
            for (int x=0;x<3;++x)
                for (int p=0;p<dk;++p){ const size_t gk=c.shells[sk].lc2pg[p];
                    for (int q=0;q<dl;++q){ const size_t gl=c.shells[sl].lc2pg[q];
                        acc[x]+=A(gk,gl)*c.scale[gk]*c.scale[gl]*buf[p+dk*(q+dl*x)]; } }
            g[c.shells[sk].atom]-=2.0*rvec3_t(acc[0], acc[1], acc[2]);   // ordered (k, l): the k derivative; the l one is the same sum relabelled (A symmetric)
        }
    return g;
}

// --- the closed-form PP gradients -----------------------------------------------------------------------
// A species field's Gaussian as an AUXILIARY shell: one primitive (coefficient 1, i.e. the raw e^{-α r²}) of momentum L at atom gAtom,
// with the per-component weights the CartGaussian's terms carry, in libcint's Cartesian order.  libcint folds CINTcommon_fac_sp(L) into
// every integral over an s or p shell, so the weights take it back out.
namespace
{
struct FieldShell { ShellDesc shell; std::vector<double> w; };
FieldShell MakeFieldShell(const Math::CartGaussian& g, size_t gAtom, const std::vector<AtomRec>& atoms)
{
    if (gAtom>=atoms.size()) throw std::out_of_range("PG_LibCint: the field's atom is past the structure");
    const AtomRec& a=atoms[gAtom];
    if (std::fabs(a.x-g.center.x)+std::fabs(a.y-g.center.y)+std::fabs(a.z-g.center.z)>1e-8)
        throw std::invalid_argument("PG_LibCint: the field Gaussian is not centred on atom "+std::to_string(gAtom));
    if (g.terms.empty()) throw std::invalid_argument("PG_LibCint: an empty field Gaussian");
    const int L=g.terms[0].p.n+g.terms[0].p.l+g.terms[0].p.m;
    FieldShell f;
    f.shell.atom=(int)gAtom; f.shell.l=L; f.shell.exp={g.alpha}; f.shell.coef={1.0};
    f.w.assign(ncart(L), 0.0);
    for (const auto& t : g.terms)
    {
        if (t.p.n+t.p.l+t.p.m!=L) throw std::invalid_argument("PG_LibCint: a field Gaussian must be homogeneous in degree (one shell)");
        int k=0;   // libcint's Cartesian order: lx descending, then ly descending
        for (int lx=L; lx>=0; --lx) for (int ly=L-lx; ly>=0; --ly) { if (lx==t.p.n && ly==t.p.l) f.w[k]+=t.c; ++k; }
    }
    for (double& x : f.w) x/=CINTcommon_fac_sp(L);
    return f;
}
} // namespace

std::vector<rvec3_t> NR_Evaluator::OverlapWithGradient(const Math::CartGaussian& g, size_t gAtom, const rvec_t& w) const
{
    const Imp& a=*itsImp;
    if (a.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the PP gradients are Cartesian-only");
    if (w.size()!=a.N) throw std::invalid_argument("PG_LibCint::OverlapWithGradient: "+std::to_string(w.size())+" weights for "+std::to_string(a.N)+" functions");
    const FieldShell f=MakeFieldShell(g, gAtom, a.atoms);
    std::vector<const ShellDesc*> sh; for (auto& s:a.shells) sh.push_back(&s);
    const int fOff=(int)sh.size(); sh.push_back(&f.shell);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(a.atoms, sh, atm,bas,env);
    const int natm=(int)a.atoms.size(), nbas=(int)sh.size(), dk=ncart(f.shell.l);
    std::vector<rvec3_t> out(a.atoms.size(), rvec3_t(0,0,0));
    std::vector<double> buf;
    for (size_t si=0;si<a.shells.size();++si)
    {
        const int di=a.dim(a.shells[si].l);
        buf.assign((size_t)di*dk*3, 0.0);
        int shls[2]={(int)si, fOff};
        int1e_ipovlp_cart(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);   // <∇i|g>
        double acc[3]={0,0,0};
        for (int x=0;x<3;++x)
            for (int p=0;p<di;++p){ const size_t gi=a.shells[si].lc2pg[p]; const double wi=w[gi]*a.scale[gi];
                for (int k=0;k<dk;++k) acc[x]+=wi*f.w[k]*buf[p+di*(k+dk*x)]; }
        const rvec3_t v(acc[0],acc[1],acc[2]);
        out[a.shells[si].atom]-=v;   // ∂φ_i/∂R_i = −∇φ_i
        out[gAtom]+=v;               // and g moves the other way (translation invariance of a 2-centre integral)
    }
    return out;
}

std::vector<rvec3_t> NR_Evaluator::OperatorThreeCGradient(const Math::CartGaussian& g, size_t gAtom, const rmat_t& D) const
{
    const Imp& a=*itsImp;
    if (a.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the PP gradients are Cartesian-only");
    if (D.rows()!=a.N || D.columns()!=a.N)
        throw std::invalid_argument("PG_LibCint::OperatorThreeCGradient: D is "+std::to_string(D.rows())+"x"+std::to_string(D.columns())+" for "+std::to_string(a.N)+" functions");
    const FieldShell f=MakeFieldShell(g, gAtom, a.atoms);
    std::vector<const ShellDesc*> sh; for (auto& s:a.shells) sh.push_back(&s);
    const int fOff=(int)sh.size(); sh.push_back(&f.shell);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(a.atoms, sh, atm,bas,env);
    const int natm=(int)a.atoms.size(), nbas=(int)sh.size(), dk=ncart(f.shell.l);
    std::vector<rvec3_t> out(a.atoms.size(), rvec3_t(0,0,0));
    std::vector<double> buf;
    for (size_t si=0;si<a.shells.size();++si)
        for (size_t sj=0;sj<a.shells.size();++sj)
        {
            const int di=a.dim(a.shells[si].l), dj=a.dim(a.shells[sj].l);
            buf.assign((size_t)di*dj*dk*3, 0.0);
            int shls[3]={(int)si,(int)sj,fOff};
            int3c1e_ip1_cart(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);   // <∇i j|g>
            double acc[3]={0,0,0};
            for (int x=0;x<3;++x)
                for (int p=0;p<di;++p){ const size_t gi=a.shells[si].lc2pg[p];
                    for (int q=0;q<dj;++q){ const size_t gj=a.shells[sj].lc2pg[q]; const double dij=D(gi,gj)*a.scale[gi]*a.scale[gj];
                        for (int k=0;k<dk;++k) acc[x]+=dij*f.w[k]*buf[p+di*(q+dj*(k+dk*x))]; } }
            const rvec3_t v=2.0*rvec3_t(acc[0],acc[1],acc[2]);   // D symmetric: the j-centre's derivative is the same sum relabelled
            out[a.shells[si].atom]-=v;
            out[gAtom]+=v;                                       // translation invariance: the field takes the opposite
        }
    return out;
}

std::vector<rvec3_t> NR_Evaluator::GaussianRepulsionGradient(size_t gAtom, double alpha, const rmat_t& D) const
{
    const Imp& a=*itsImp;
    if (a.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the PP gradients are Cartesian-only");
    if (D.rows()!=a.N || D.columns()!=a.N)
        throw std::invalid_argument("PG_LibCint::GaussianRepulsionGradient: D is "+std::to_string(D.rows())+"x"+std::to_string(D.columns())+" for "+std::to_string(a.N)+" functions");
    if (gAtom>=a.atoms.size()) throw std::out_of_range("PG_LibCint::GaussianRepulsionGradient: the atom is past the structure");
    const AtomRec& at=a.atoms[gAtom];
    const FieldShell f=MakeFieldShell(Math::CartGaussian{rvec3_t(at.x,at.y,at.z), alpha, {Math::CartTerm{Math::Monomial{0,0,0},1.0}}}, gAtom, a.atoms);
    std::vector<const ShellDesc*> sh; for (auto& s:a.shells) sh.push_back(&s);
    const int fOff=(int)sh.size(); sh.push_back(&f.shell);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(a.atoms, sh, atm,bas,env);
    const int natm=(int)a.atoms.size(), nbas=(int)sh.size();
    std::vector<rvec3_t> out(a.atoms.size(), rvec3_t(0,0,0));
    std::vector<double> buf;
    for (size_t si=0;si<a.shells.size();++si)
        for (size_t sj=0;sj<a.shells.size();++sj)
        {
            const int di=a.dim(a.shells[si].l), dj=a.dim(a.shells[sj].l);
            int shls[3]={(int)si,(int)sj,fOff};
            for (int which=0; which<2; ++which)
            {
                buf.assign((size_t)di*dj*3, 0.0);
                (which==0 ? int3c2e_ip1_cart : int3c2e_ip2_cart)(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
                double acc[3]={0,0,0};
                for (int x=0;x<3;++x)
                    for (int p=0;p<di;++p){ const size_t gi=a.shells[si].lc2pg[p];
                        for (int q=0;q<dj;++q){ const size_t gj=a.shells[sj].lc2pg[q];
                            acc[x]+=D(gi,gj)*a.scale[gi]*a.scale[gj]*f.w[0]*buf[p+di*(q+dj*x)]; } }
                const int on = which==0 ? a.shells[si].atom : (int)gAtom;
                out[on]-=(which==0 ? 2.0 : 1.0)*rvec3_t(acc[0],acc[1],acc[2]);   // ip1 on i (and, D symmetric, on j); ip2 on the charge
            }
        }
    return out;
}

// THE HF GRADIENT.  f = Σ Da_ij Db_kl (ij|kl) (Coulomb) or (ik|jl) (exchange), i, j on this evaluator, k, l on the other.  int2e_ip1
// differentiates the FIRST function; the symmetries (ij|kl) = (ji|kl) = (kl|ij) and real orbitals give every other position from a
// relabelling, with Da, Db symmetric:
//   Coulomb : the i and j derivatives are equal, Bra(a; Da, Db); the k and l ones are Bra(b; Db, Da);
//   exchange: (ik|jl) = (jl|ik) makes j equal i, BraX(a; Da, Db); (ik|jl) = (ki|lj) makes k (and l) BraX(b; Db, Da),
// each counted twice.  ∂φ/∂R = −∇φ.
std::vector<rvec3_t> NR_Evaluator::ERI4Gradient(const NR_Evaluator& other, const rmat_t& Da, const rmat_t& Db, bool exchange) const
{
    const Imp& A=*itsImp; const Imp& B=*other.itsImp;
    if (A.spherical || B.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the ERI4 gradient is Cartesian-only");
    if (Da.rows()!=A.N || Da.columns()!=A.N || Db.rows()!=B.N || Db.columns()!=B.N)
        throw std::invalid_argument("PG_LibCint::ERI4Gradient: Da is "+std::to_string(Da.rows())+" and Db "+std::to_string(Db.rows())+" square, for "
                                    +std::to_string(A.N)+" and "+std::to_string(B.N)+" functions");
    if (A.atoms.size()!=B.atoms.size()) throw std::invalid_argument("PG_LibCint::ERI4Gradient: the two blocks are on different structures");
    std::vector<const ShellDesc*> sh; for (auto& s:A.shells) sh.push_back(&s);
    const size_t oB=sh.size();        for (auto& s:B.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(A.atoms, sh, atm,bas,env);
    const int natm=(int)A.atoms.size(), nbas=(int)sh.size();
    std::vector<rvec3_t> g(A.atoms.size(), rvec3_t(0,0,0));
    std::vector<double> buf;
    // One side: X's shells carry the differentiated function (offset oX), Y's the other density (offset oY).
    auto side=[&](const Imp& X, size_t oX, const rmat_t& DX, const Imp& Y, size_t oY, const rmat_t& DY)
    {
        for (size_t si=0;si<X.shells.size();++si)
        {
            double acc[3]={0,0,0};
            const int di=X.dim(X.shells[si].l);
            for (size_t sj=0;sj<X.shells.size();++sj)
                for (size_t sk=0;sk<Y.shells.size();++sk)
                    for (size_t sl=0;sl<Y.shells.size();++sl)
                    {
                        const int dj=X.dim(X.shells[sj].l), dk=Y.dim(Y.shells[sk].l), dl=Y.dim(Y.shells[sl].l);
                        buf.assign((size_t)di*dj*dk*dl*3, 0.0);
                        // Coulomb (i j | k l); exchange (i k | j l): the buffer's index order follows the shell order.
                        int shls[4];
                        if (!exchange) { shls[0]=int(oX+si); shls[1]=int(oX+sj); shls[2]=int(oY+sk); shls[3]=int(oY+sl); }
                        else           { shls[0]=int(oX+si); shls[1]=int(oY+sk); shls[2]=int(oX+sj); shls[3]=int(oY+sl); }
                        int2e_ip1_cart(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
                        for (int x=0;x<3;++x)
                            for (int p=0;p<di;++p){ const size_t gi=X.shells[si].lc2pg[p];
                                for (int q=0;q<dj;++q){ const size_t gj=X.shells[sj].lc2pg[q]; const double dij=DX(gi,gj)*X.scale[gi]*X.scale[gj];
                                    if (dij==0.0) continue;
                                    for (int r=0;r<dk;++r){ const size_t gk=Y.shells[sk].lc2pg[r];
                                        for (int t=0;t<dl;++t){ const size_t gl=Y.shells[sl].lc2pg[t];
                                            const size_t at = !exchange ? p+di*(q+dj*(r+dk*(t+dl*x))) : p+di*(r+dk*(q+dj*(t+dl*x)));
                                            acc[x]+=dij*DY(gk,gl)*Y.scale[gk]*Y.scale[gl]*buf[at]; } } } }
                    }
            g[X.shells[si].atom]-=2.0*rvec3_t(acc[0],acc[1],acc[2]);   // twice: the partner position on the same side
        }
    };
    side(A, 0,  Da, B, oB, Db);
    side(B, oB, Db, A, 0,  Da);
    return g;
}

rsmat_t NR_Evaluator::RepulsionTwoC_Matrix() const
{
    const Imp& c=*itsImp;
    if (c.spherical) throw std::logic_error("PG_LibCint::NR_Evaluator: the 2-centre Coulomb metric is Cartesian-only");
    std::vector<const ShellDesc*> sh; for (auto& s:c.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(c.atoms, sh, atm,bas,env);
    const int natm=(int)c.atoms.size(), nbas=(int)sh.size();
    rsmat_t M(c.N);
    std::vector<double> buf;
    for (size_t sk=0;sk<c.shells.size();++sk)
        for (size_t sl=sk;sl<c.shells.size();++sl)
        {
            const int dk=c.dim(c.shells[sk].l), dl=c.dim(c.shells[sl].l);
            buf.assign((size_t)dk*dl, 0.0);
            int shls[2]={(int)sk,(int)sl};
            int2c2e_cart(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
            for (int p=0;p<dk;++p){ const size_t gk=c.shells[sk].lc2pg[p];
                for (int q=0;q<dl;++q){ const size_t gl=c.shells[sl].lc2pg[q];
                    M(gk,gl)=buf[p+dk*q]*c.scale[gk]*c.scale[gl]; } }
        }
    return M;
}

Projector3<double> NR_Evaluator::OverlapThreeC_Matrix  (const NR_Evaluator& fit) const {return Build3(fit,0);}
Projector3<double> NR_Evaluator::RepulsionThreeC_Matrix(const NR_Evaluator& fit) const {return Build3(fit,1);}

// --- 4-centre (ab|cd): dense reordered+renormalized tensor over four (this/partner) slots ----------
std::vector<double> NR_Evaluator::Compute4(const NR_Evaluator& B, const NR_Evaluator& C,
                                           const NR_Evaluator& D) const
{
    const Imp& A=*itsImp; const Imp& mb=*B.itsImp; const Imp& mc=*C.itsImp; const Imp& md=*D.itsImp;
    size_t Na=A.N,Nb=mb.N,Nc=mc.N,Nd=md.N;
    std::vector<const ShellDesc*> sh;
    size_t oA=0;            for (auto& s:A.shells ) sh.push_back(&s);
    size_t oB=sh.size();    for (auto& s:mb.shells) sh.push_back(&s);
    size_t oC=sh.size();    for (auto& s:mc.shells) sh.push_back(&s);
    size_t oD=sh.size();    for (auto& s:md.shells) sh.push_back(&s);
    std::vector<int> atm,bas; std::vector<double> env; Imp::BuildMol(A.atoms, sh, atm,bas,env);
    int natm=(int)A.atoms.size(), nbas=(int)sh.size();

    std::vector<double> G(Na*Nb*Nc*Nd, 0.0);
    auto idx=[&](size_t ia,size_t ib,size_t ic,size_t id){return ((ia*Nb+ib)*Nc+ic)*Nd+id;};
    std::vector<double> buf;
    for (size_t sa=0;sa<A.shells.size();++sa)
     for (size_t sb=0;sb<mb.shells.size();++sb)
      for (size_t sc=0;sc<mc.shells.size();++sc)
       for (size_t sd=0;sd<md.shells.size();++sd)
       {
           int da=A.dim(A.shells[sa].l), db=mb.dim(mb.shells[sb].l),
               dc=mc.dim(mc.shells[sc].l), dd=md.dim(md.shells[sd].l);
           buf.assign((size_t)da*db*dc*dd, 0.0);
           int shls[4]={(int)(oA+sa),(int)(oB+sb),(int)(oC+sc),(int)(oD+sd)};
           IntFn fn = A.spherical ? int2e_sph : int2e_cart;
           fn(buf.data(), nullptr, shls, atm.data(), natm, bas.data(), nbas, env.data(), nullptr, nullptr);
           for (int p=0;p<da;++p){ size_t gia=A.shells[sa].lc2pg[p];
            for (int q=0;q<db;++q){ size_t gib=mb.shells[sb].lc2pg[q];
             for (int rr=0;rr<dc;++rr){ size_t gic=mc.shells[sc].lc2pg[rr];
              for (int s4=0;s4<dd;++s4){ size_t gid=md.shells[sd].lc2pg[s4];
                G[idx(gia,gib,gic,gid)] = buf[p+da*(q+db*(rr+dc*s4))]
                    * A.scale[gia]*mb.scale[gib]*mc.scale[gic]*md.scale[gid]; } } } }
       }
    return G;
}

ERI4 NR_Evaluator::DirectMatrix(const NR_Evaluator& p) const
{
    size_t Na=itsImp->N, Nc=p.itsImp->N;
    std::vector<double> G=Compute4(*this, p, p);   // (a a | c c): A=this,B=this,C=p,D=p
    auto g=[&](size_t ia,size_t ib,size_t ic,size_t id){return G[((ia*Na+ib)*Nc+ic)*Nc+id];};
    ERI4 J(Na,Nc);
    for (size_t ia=0;ia<Na;++ia) for (size_t ib=ia;ib<Na;++ib)
    {
        rsmat_t& Jab=J(ia,ib);
        for (size_t ic=0;ic<Nc;++ic) for (size_t id=ic;id<Nc;++id) Jab(ic,id)=g(ia,ib,ic,id);
    }
    return J;
}

ERI4 NR_Evaluator::ExchangeMatrix(const NR_Evaluator& pb) const
{
    size_t Na=itsImp->N, Nb=pb.itsImp->N;
    std::vector<double> G=Compute4(pb, *this, pb);   // slots (a b | a b): A=this,B=pb,C=this,D=pb
    auto g=[&](size_t ia,size_t ib,size_t ic,size_t id){return G[((ia*Nb+ib)*Na+ic)*Nb+id];};
    ERI4 K(Na,Nb);
    for (size_t ia=0;ia<Na;++ia) for (size_t ib=0;ib<Nb;++ib) for (size_t ic=ia;ic<Na;++ic)
    {
        rsmat_t& Kac=K(ia,ic);
        for (size_t id=0;id<Nb;++id)
        {
            double v=g(ia,ib,ic,id);
            if (ib==id)     Kac(ib,id) = v;
            else if (ib<id) Kac(ib,id)+= 0.5*v;
            else            Kac(id,ib)+= 0.5*v;
        }
    }
    return K;
}

} //namespace

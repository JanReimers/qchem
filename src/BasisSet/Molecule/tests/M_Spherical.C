// File: BasisSet/Molecule/tests/M_Spherical.C
//
// Increment-1 guard for PG_Spherical_MnD: the real-solid-harmonic transform must produce mutually
// orthogonal spherical functions.  We assemble the raw spherical self-overlap <chi_m|chi_m'> from the
// trusted Cartesian overlap (GaussianRF::Overlap2C) and the transform coefficients, and require it to be
// diagonal -- which is exactly the statement that the 2l+1 functions of a shell are orthogonal.  This
// validates the transform end to end against already-oracle-verified Cartesian integrals.
#include "gtest/gtest.h"
#include <algorithm>   // std::sort (the ordering-invariant spectra in the view gates)
#include <cmath>
#include <complex>     // std::real (the pencil spectra)
#include <vector>

import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;      // GaussianRF (Cartesian kernels)
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;    // Polarization
import qchem.Math.Angular;                                            // SphericalShell, CartTerm
import qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD;            // NR_Evaluator (1E)
import qchem.BasisSet;                                                // Real_BS / Real_OIBS iteration (view gates)
import qchem.BasisSet.Orbital_1E_IBS;                                 // Real_OIBS faces
import qchem.BasisSet.Molecule.PG_Cart;                               // the wrapped Cartesian family
import qchem.BasisSet.Molecule.PG_Spherical;                          // the native spherical family (oracle)
import qchem.BasisSet.Molecule.PG_Spherical.LatticeView;              // MakeSphericalLatticeView (under test)
import qchem.Structure;                                               // Molecule, Atom
import qchem.Blaze;                                                   // blazem::eigen/inv/dot (view gates)
import qchem.Types;
using namespace qchem;

using namespace qchem::BasisSet::Molecule::Evaluators;
using qchem::Math::SphericalShell;
using PG_Spherical_MnD::NR_Evaluator;
using PG_Cart_MnD::GaussianRF;

// Raw <chi_i | chi_j> over a single-primitive shell of total angular momentum l, centred at the origin.
static void check_orthogonal(int l)
{
    GaussianRF rf(1.3, rvec3_t(0,0,0), l);     // one exponent
    auto sh = SphericalShell(l);
    ASSERT_EQ(sh.size(), size_t(2*l+1));

    for (size_t i=0;i<sh.size();++i)
        for (size_t j=0;j<sh.size();++j)
        {
            double s = 0.0;
            for (const auto& ta : sh[i])
                for (const auto& tb : sh[j])
                    s += ta.c * tb.c * rf.Overlap2C(rf, ta.p, tb.p);
            if (i!=j) EXPECT_NEAR(s, 0.0, 1e-12) << "l="<<l<<" components ("<<i<<","<<j<<") not orthogonal";
            else      EXPECT_GT(std::fabs(s), 1e-12) << "l="<<l<<" component "<<i<<" has zero norm";
        }
}

TEST(M_Spherical, p_orthogonal) { check_orthogonal(1); }
TEST(M_Spherical, d_orthogonal) { check_orthogonal(2); }
TEST(M_Spherical, f_orthogonal) { check_orthogonal(3); }

// --- Increment 2: the 1E spherical evaluator ------------------------------------------------------
// Two centres on the z-axis, each carrying one s shell and one d shell.  Building the evaluator's
// Overlap/Grad2 through the transform-on-Cartesian kernels must give: a normalised diagonal (=1), an
// orthonormal block per centre (s vs d on one centre -> 0, the 5 d's mutually -> 0), and symmetry.
TEST(M_Spherical, evaluator_1E)
{
    std::vector<GaussianRF> radials;          // own the radials; reserve so the pointers stay valid
    radials.reserve(4);
    radials.push_back(GaussianRF(0.50, rvec3_t(0,0, 0.0), 0));   // s on centre A
    radials.push_back(GaussianRF(0.80, rvec3_t(0,0, 0.0), 2));   // d on centre A
    radials.push_back(GaussianRF(0.50, rvec3_t(0,0, 1.4), 0));   // s on centre B
    radials.push_back(GaussianRF(0.80, rvec3_t(0,0, 1.4), 2));   // d on centre B

    NR_Evaluator ev;
    std::vector<int> centre;                  // which centre each component belongs to (0=A, 1=B)
    auto addShell = [&](size_t r, int l, int c)
    {
        for (auto& terms : SphericalShell(l)) { ev.comps.push_back({&radials[r], terms}); centre.push_back(c); }
    };
    addShell(0,0,0); addShell(1,2,0);   // A: 1 + 5 = 6 components
    addShell(2,0,1); addShell(3,2,1);   // B: 1 + 5 = 6 components
    ev.Init();
    ASSERT_EQ(ev.size(), size_t(12));

    for (size_t i=0;i<ev.size();++i)
        for (size_t j=0;j<ev.size();++j)
        {
            double s = ev.Overlap(i,j);
            EXPECT_NEAR(s, ev.Overlap(j,i), 1e-13) << "overlap not symmetric ("<<i<<","<<j<<")";
            EXPECT_NEAR(ev.Grad2(i,j), ev.Grad2(j,i), 1e-12) << "grad2 not symmetric ("<<i<<","<<j<<")";
            if (i==j)
                EXPECT_NEAR(s, 1.0, 1e-12) << "component "<<i<<" not normalised";
            else if (centre[i]==centre[j])
                EXPECT_NEAR(s, 0.0, 1e-12) << "same-centre components ("<<i<<","<<j<<") not orthogonal";
        }
}

// --- Increment 3: the 3-centre (DFT) and 4-centre (HF) kernels ------------------------------------
// Two d shells on two centres (10 components).  The transform-summed 3C/4C must obey the standard
// integral permutation symmetries -- (ab|c)=(ba|c) and (ab|cd)=(ba|cd)=(ab|dc)=(cd|ab) -- which the
// Cartesian kernels (oracle-verified) already satisfy and the spherical transform must preserve.
TEST(M_Spherical, evaluator_3C4C_symmetry)
{
    std::vector<GaussianRF> radials; radials.reserve(2);
    radials.push_back(GaussianRF(0.70, rvec3_t(0,0, 0.0), 2));   // d on A
    radials.push_back(GaussianRF(0.60, rvec3_t(0,0, 1.2), 2));   // d on B

    NR_Evaluator ev;
    auto addShell = [&](size_t r,int l){ for (auto& t : SphericalShell(l)) ev.comps.push_back({&radials[r],t}); };
    addShell(0,2); addShell(1,2);     // 10 components (0..4 on A, 5..9 on B)
    ev.Init();
    ASSERT_EQ(ev.size(), size_t(10));

    const int idx[][4] = { {0,1,2,3}, {0,5,7,9}, {2,2,8,8}, {1,6,3,4} };
    for (const auto& q : idx)
    {
        int a=q[0], b=q[1], c=q[2], d=q[3];
        double abcd = ev.FourC(a,ev,b,ev,c,ev,d);
        EXPECT_NEAR(abcd, ev.FourC(b,ev,a,ev,c,ev,d), 1e-11) << "(ab|cd) != (ba|cd)";
        EXPECT_NEAR(abcd, ev.FourC(a,ev,b,ev,d,ev,c), 1e-11) << "(ab|cd) != (ab|dc)";
        EXPECT_NEAR(abcd, ev.FourC(c,ev,d,ev,a,ev,b), 1e-11) << "(ab|cd) != (cd|ab)";

        EXPECT_NEAR(ev.OverlapThreeC  (a,ev,b,ev,c), ev.OverlapThreeC  (b,ev,a,ev,c), 1e-12) << "overlap3c (ab|c)!=(ba|c)";
        EXPECT_NEAR(ev.RepulsionThreeC(a,ev,b,ev,c), ev.RepulsionThreeC(b,ev,a,ev,c), 1e-12) << "repulsion3c (ab|c)!=(ba|c)";
    }
}

// --- Increment 5: the DFT fit-basis kernels (Charge + 2-centre Coulomb metric) --------------------
// The defining property of a spherical fit basis: every l>0 real solid harmonic is contaminant-free, so
// its monopole charge integral is exactly zero -- only the s components carry charge.  (This is what makes
// the spherical fit smaller than the Cartesian one: the Cartesian d-shell's xx+yy+zz is an l=0 charge-
// carrying function the spherical shell drops.)  Plus the Coulomb metric is symmetric + positive.
TEST(M_Spherical, fit_kernels)
{
    std::vector<GaussianRF> radials; radials.reserve(2);
    radials.push_back(GaussianRF(0.90, rvec3_t(0,0,0), 0));   // s
    radials.push_back(GaussianRF(0.70, rvec3_t(0,0,0), 2));   // d

    NR_Evaluator ev;
    auto add = [&](size_t r,int l){ for (auto& t : SphericalShell(l)) ev.comps.push_back({&radials[r],t}); };
    add(0,0); add(1,2);     // 1 s + 5 d = 6 components
    ev.Init();

    EXPECT_GT(ev.Charge(0), 0.0) << "s component should carry charge";
    for (size_t i=1;i<6;i++)
        EXPECT_NEAR(ev.Charge(i), 0.0, 1e-12) << "d harmonic " << i << " is not contaminant-free (charge!=0)";

    for (size_t i=0;i<6;i++)
    {
        EXPECT_GT(ev.Repulsion2C(i,i), 0.0) << "self Coulomb metric must be positive (" << i << ")";
        for (size_t j=0;j<6;j++)
            EXPECT_NEAR(ev.Repulsion2C(i,j), ev.Repulsion2C(j,i), 1e-12) << "Coulomb metric not symmetric ("<<i<<","<<j<<")";
    }
}

// --- The SPHERICAL LATTICE VIEW (doc/SphericalLatticePlan.md I1) --------------------------------
// The view wraps a Cartesian basis in its contaminant-free spherical span through the abstract
// faces alone (engine-blind).  Two gates:
//  (1) NULL TEST: an s/p-only basis has NO contaminants, so the view must reproduce the inner
//      matrices EXACTLY (T = identity).
//  (2) NATIVE-FAMILY ORACLE: a d shell viewed spherically must match the NATIVE PG_Spherical
//      family (same exponents) matrix-for-matrix -- the m-ordering and normalisation conventions
//      agree end to end, or this fails loudly.
TEST(M_SphericalView, sp_null_test)
{
    Molecule mol;
    mol.Insert(new Atom(8, 0, Vector3D<double>(0,0,0)));
    mol.Insert(new Atom(8, 0, Vector3D<double>(0,0,1.8)));
    auto* cart=new qchem::BasisSet::Molecule::PG_Cart::BasisSet;
    cart->Insert(new qchem::BasisSet::Molecule::PG_Cart::Orbital_IBS(rvec_t{0.5,1.7}, 1, &mol));
    std::shared_ptr<const qchem::BasisSet::Real_BS> cbs(cart);
    auto view=qchem::BasisSet::Molecule::PG_Spherical::MakeSphericalLatticeView(cbs);

    const qchem::BasisSet::Real_OIBS* vi=nullptr;
    for (auto b : const_cast<qchem::BasisSet::Real_BS&>(*view).Iterate<qchem::BasisSet::Real_OIBS>()) { vi=b; break; }
    const qchem::BasisSet::Real_OIBS* ci=nullptr;
    for (auto b : const_cast<qchem::BasisSet::Real_BS&>(*cbs).Iterate<qchem::BasisSet::Real_OIBS>()) { ci=b; break; }
    ASSERT_TRUE(vi && ci);
    ASSERT_EQ(vi->GetNumFunctions(), ci->GetNumFunctions());       // no contaminants to drop
    const auto Sv=vi->MakeOverlap(), Sc=ci->MakeOverlap();
    const auto Tv=vi->MakeKinetic(), Tc2=ci->MakeKinetic();
    const auto Vv=vi->MakeNuclear(&mol), Vc=ci->MakeNuclear(&mol);
    const size_t n=ci->GetNumFunctions();
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            EXPECT_NEAR(Sv(i,j), Sc(i,j), 1e-13);
            EXPECT_NEAR(Tv(i,j), Tc2(i,j), 1e-12);
            EXPECT_NEAR(Vv(i,j), Vc(i,j), 1e-12);
        }
}

TEST(M_SphericalView, d_matches_native_family)
{
    Molecule mol;
    mol.Insert(new Atom(25, 0, Vector3D<double>(0,0,0)));
    mol.Insert(new Atom(25, 0, Vector3D<double>(0,0,2.1)));
    const rvec_t ex{0.8,2.4};

    auto* cart=new qchem::BasisSet::Molecule::PG_Cart::BasisSet;
    cart->Insert(new qchem::BasisSet::Molecule::PG_Cart::Orbital_IBS(ex, 2, &mol));
    std::shared_ptr<const qchem::BasisSet::Real_BS> cbs(cart);
    auto view=qchem::BasisSet::Molecule::PG_Spherical::MakeSphericalLatticeView(cbs);

    auto* nat=new qchem::BasisSet::Molecule::PG_Spherical::BasisSet;
    nat->Insert(new qchem::BasisSet::Molecule::PG_Spherical::Orbital_IBS(ex, 2, &mol));
    std::shared_ptr<const qchem::BasisSet::Real_BS> nbs(nat);

    const qchem::BasisSet::Real_OIBS* vi=nullptr;
    for (auto b : const_cast<qchem::BasisSet::Real_BS&>(*view).Iterate<qchem::BasisSet::Real_OIBS>()) { vi=b; break; }
    const qchem::BasisSet::Real_OIBS* ni=nullptr;
    for (auto b : const_cast<qchem::BasisSet::Real_BS&>(*nbs).Iterate<qchem::BasisSet::Real_OIBS>()) { ni=b; break; }
    ASSERT_TRUE(vi && ni);
    ASSERT_EQ(vi->GetNumFunctions(), ni->GetNumFunctions());       // contaminants dropped == native count
    // The two families ORDER their functions differently (view = the inner Cartesian layout with each
    // d 6-block compressed to 5; native = its own layout), so matrices are compared through
    // ordering-INVARIANT probes of the same span:
    //   (i)  sorted eigenvalues of S;
    //   (ii) the pencil spectra eig(S^{-1}T), eig(S^{-1}V) -- basis-change invariants of the operators;
    //   (iii) the span projector f(r)^T S^{-1} f(r) at sample points (invariant under ANY in-span
    //         basis change -- the pointwise "is it the same function space" probe).
    const size_t n=ni->GetNumFunctions();
    const rmat_t Sv(vi->MakeOverlap()), Sn(ni->MakeOverlap());
    const rmat_t Tv(vi->MakeKinetic()), Tn(ni->MakeKinetic());
    const rmat_t Vv(vi->MakeNuclear(&mol)), Vn(ni->MakeNuclear(&mol));
    auto sortedEigsSym=[n](rmat_t A)
    {
        hmat_t<double> H(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) H(i,j)=0.5*(A(i,j)+A(j,i));
        rvec_t w; blazem::eigen(H,w);
        std::vector<double> s(w.size());                        // explicit copy: std::sort + blaze
        for (size_t i=0;i<w.size();i++) s[i]=w[i];              // iterators do not mix (CLAUDE.md)
        std::sort(s.begin(),s.end());
        return s;
    };
    auto pencilEigs=[n](const rmat_t& A, const rmat_t& S)
    {
        rmat_t M=blazem::inv(S)*A;                    // similar to S^{-1/2} A S^{-1/2}: real spectrum
        cvec_t w; blazem::eigen(M,w);
        std::vector<double> s(w.size());
        for (size_t i=0;i<w.size();i++) s[i]=std::real(w[i]);
        std::sort(s.begin(),s.end());
        return s;
    };
    const auto es_v=sortedEigsSym(Sv), es_n=sortedEigsSym(Sn);
    const auto et_v=pencilEigs(Tv,Sv), et_n=pencilEigs(Tn,Sn);
    const auto ev_v=pencilEigs(Vv,Sv), ev_n=pencilEigs(Vn,Sn);
    for (size_t i=0;i<n;i++)
    {
        EXPECT_NEAR(es_v[i], es_n[i], 1e-11) << "eig(S)["<<i<<"]";
        EXPECT_NEAR(et_v[i], et_n[i], 1e-9)  << "eig(S^-1 T)["<<i<<"]";
        EXPECT_NEAR(ev_v[i], ev_n[i], 1e-9)  << "eig(S^-1 V)["<<i<<"]";
    }
    const rmat_t Siv=blazem::inv(Sv), Sin=blazem::inv(Sn);
    for (const rvec3_t r : {rvec3_t(0.3,0.1,0.2), rvec3_t(-0.2,0.4,1.9)})
    {
        const rvec_t fv=(*vi)(r), fn=(*ni)(r);
        const double pv=blazem::dot(fv,Siv*fv), pn=blazem::dot(fn,Sin*fn);
        EXPECT_NEAR(pv, pn, 1e-10) << "span projector at r";
    }
}

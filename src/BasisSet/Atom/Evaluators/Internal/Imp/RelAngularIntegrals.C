// File: BasisSet/Atom/Evaluators/Internal/Imp/RelAngularIntegrals.C
module;
#include <cassert>
#include <cmath>            // std::abs (double)
#include <initializer_list> // the {-0.5, 0.5} spin loops

module qchem.BasisSet.Atom.Evaluators.Internal.RelAngularIntegrals;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;
import qchem.BasisSet.Atom.Evaluators.Internal.Wigner3j; //Wigner::wigner3j / clebschGordan (home-grown)
import qchem.Symmetry.Spherical; //SphericalSpinor::j, l
import qchem.Math; //FourPi2
import qchem.Blaze;

namespace RelAngularIntegrals
{

// CG coefficient <la, mj-ms; ½, ms | ja, mj>
static double CG(int la, double ja, double mj, double ms)
{
    double mla = mj - ms;
    if (std::abs(mla) > la) return 0.0;
    return Wigner::clebschGordan(la, 0.5, ja, mla, ms, mj);
}

rvec11_t Direct (int κa, int κc, double mja, double mjc)
{
    int    la=Symmetry::SphericalSpinor::l(κa), lc=Symmetry::SphericalSpinor::l(κc);
    double ja=Symmetry::SphericalSpinor::j(κa), jc=Symmetry::SphericalSpinor::j(κc);
    assert(std::abs(mja)<=ja);
    assert(std::abs(mjc)<=jc);
    rvec11_t Ak(0.0);
    for (double ms_a : {-0.5, 0.5})
    {
        double cga=CG(la,ja,mja,ms_a);
        if (cga==0.0) continue;
        double cga2=cga*cga;
        for (double ms_c : {-0.5, 0.5})
        {
            double cgc=CG(lc,jc,mjc,ms_c);
            if (cgc==0.0) continue;
            int mla=(int)(mja-ms_a), mlc=(int)(mjc-ms_c);
            Ak += cga2 * cgc*cgc * AngularIntegrals::Direct (la,lc,mla,mlc);
        }
    }
    return Ak;
}

rvec11_t Exchange(int κa, int κb, double mja, double mjb)
{
    int    la=Symmetry::SphericalSpinor::l(κa), lb=Symmetry::SphericalSpinor::l(κb);
    double ja=Symmetry::SphericalSpinor::j(κa), jb=Symmetry::SphericalSpinor::j(κb);
    assert(std::abs(mja)<=ja);
    assert(std::abs(mjb)<=jb);
    rvec11_t Ak(0.0);
    for (double ms : {-0.5, 0.5}) //spin must match across exchange
    {
        double cga=CG(la,ja,mja,ms);
        if (cga==0.0) continue;
        double cgb=CG(lb,jb,mjb,ms);
        if (cgb==0.0) continue;
        int mla=(int)(mja-ms), mlb=(int)(mjb-ms);
        Ak += cga*cga * cgb*cgb * AngularIntegrals::Exchange(la,lb,mla,mlb);
    }
    return Ak;
}

rvec11_t Direct (int κa, int κc)
{
    double ja=Symmetry::SphericalSpinor::j(κa), jc=Symmetry::SphericalSpinor::j(κc);
    rvec11_t Ak(0.0);
    for (double mja=-ja; mja<=ja; mja+=1.0)
    for (double mjc=-jc; mjc<=jc; mjc+=1.0)
        Ak += Direct (κa,κc,mja,mjc);
    return Ak;
}

// Relativistic exchange Slater-integral coefficients for closed subshells (κa with κc).
//   Ak^k = 4π² · 2 · (ja k jc; ½ 0 -½)²   for k in [|la-lc|, la+lc] step 2.
// Using the jj 3j directly (rather than decomposing the large component into l,ml and
// CG-coupling) gives the correct k-selection: it forbids k=2 for a j=½ pair, and makes
// the s<->p cross term κ-independent and equal to the nonrelativistic value.  The factor
// 2 (not 2jc+1) keeps the coefficient symmetric in a<->c and reproduces the NR limit:
// it gives 4π² for s-s and matches DHF references for He, Be, Ne and B 2p to <0.5%.
// The earlier CG-decomposition was correct only for s-s (l=0); for l>0 it produced
// spurious multipoles and a factor-of-2-too-small cross exchange.
rvec11_t Exchange(int κa, int κc)
{
    int    la=Symmetry::SphericalSpinor::l(κa), lc=Symmetry::SphericalSpinor::l(κc);
    double ja=Symmetry::SphericalSpinor::j(κa), jc=Symmetry::SphericalSpinor::j(κc);
    rvec11_t Ak(0.0);
    for (int k=std::abs(la-lc); k<=la+lc; k+=2)
    {
        double w=Wigner::wigner3j(ja,(double)k,jc,0.5,0.0,-0.5);
        Ak[k]=FourPi2*2.0*w*w;
    }
    return Ak;
}

} // namespace

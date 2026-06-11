// File: BasisSet1/Atom/Evaluators/Internal/Imp/RelAngularIntegrals.C
module;
#include "wignerSymbols/wignerSymbols-cpp.h"
#include <blaze/Math.h>
#include <cassert>
#include <cmath>

module qchem.BasisSet.Atom.Evaluators.Internal.RelAngularIntegrals;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;
import qchem.Symmetry.Spherical; //SphericalSpinor::j, l

namespace RelAngularIntegrals
{

// CG coefficient <la, mj-ms; ½, ms | ja, mj>
static double CG(int la, double ja, double mj, double ms)
{
    double mla = mj - ms;
    if (std::abs(mla) > la) return 0.0;
    return WignerSymbols::clebschGordan(la, 0.5, ja, mla, ms, mj);
}

rvec11_t Coulomb(int κa, int κc, double mja, double mjc)
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
            Ak += cga2 * cgc*cgc * AngularIntegrals::Coulomb(la,lc,mla,mlc);
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

rvec11_t Coulomb(int κa, int κc)
{
    double ja=Symmetry::SphericalSpinor::j(κa), jc=Symmetry::SphericalSpinor::j(κc);
    rvec11_t Ak(0.0);
    for (double mja=-ja; mja<=ja; mja+=1.0)
    for (double mjc=-jc; mjc<=jc; mjc+=1.0)
        Ak += Coulomb(κa,κc,mja,mjc);
    return Ak;
}

rvec11_t Exchange(int κa, int κb)
{
    double ja=Symmetry::SphericalSpinor::j(κa), jb=Symmetry::SphericalSpinor::j(κb);
    rvec11_t Ak(0.0);
    for (double mja=-ja; mja<=ja; mja+=1.0)
    for (double mjb=-jb; mjb<=jb; mjb+=1.0)
        Ak += Exchange(κa,κb,mja,mjb);
    return Ak;
}

} // namespace

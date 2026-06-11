// File: BasisSet1/Atom/Evaluators/Internal/Imp/RelWigner3j.C
module;
#include "wignerSymbols/wignerSymbols-cpp.h"
#include <cassert>
#include <iostream>
module qchem.BasisSet.Atom.Evaluators.Internal.RelWigner3j;
import qchem.Symmetry.Spherical; //For SphericalSpinor::j(κ)

RelWigner3j RelWigner3j::w3j;

RelWigner3j::RelWigner3j()
{
    std::cout << "Initializing relativistic Wigner 3j tables LMax=" << LMax << std::endl;

    for (int κa=-(KMax); κa<=KMax; κa++)
    {
        if (κa==0) continue;
        double ja=Symmetry::SphericalSpinor::j(κa);
        for (int κb=-(KMax); κb<=KMax; κb++)
        {
            if (κb==0) continue;
            double jb=Symmetry::SphericalSpinor::j(κb);
            for (int k=0; k<=KkMax; k++)
            {
                Data[κ_idx(κa)][κ_idx(κb)][k] =
                    WignerSymbols::wigner3j(ja, k, jb, 0.5, 0.0, -0.5);

                for (double mja=-ja; mja<=ja; mja+=1.0)
                    for (double mjb=-jb; mjb<=jb; mjb+=1.0)
                        Data_m[κ_idx(κa)][κ_idx(κb)][k][mj_idx(mja)][mj_idx(mjb)] =
                            WignerSymbols::wigner3j(ja, k, jb, -mja, mja-mjb, mjb);
            }
        }
    }
}

double RelWigner3j::operator()(int κa, int κb, int k) const
{
    assert(κa!=0 && abs(κa)<=KMax);
    assert(κb!=0 && abs(κb)<=KMax);
    assert(k>=0 && k<=KkMax);
    return Data[κ_idx(κa)][κ_idx(κb)][k];
}

double RelWigner3j::operator()(int κa, int κb, int k, double mja, double mjb) const
{
    assert(κa!=0 && abs(κa)<=KMax);
    assert(κb!=0 && abs(κb)<=KMax);
    assert(k>=0 && k<=KkMax);
    double ja=Symmetry::SphericalSpinor::j(κa);
    double jb=Symmetry::SphericalSpinor::j(κb);
    assert(mja>=-ja && mja<=ja);
    assert(mjb>=-jb && mjb<=jb);
    return Data_m[κ_idx(κa)][κ_idx(κb)][k][mj_idx(mja)][mj_idx(mjb)];
}

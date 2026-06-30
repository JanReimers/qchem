// File: BasisSet/Molecule/Evaluators/PG_Spherical_MnD/Imp/SolidHarmonics.C
//
// Real regular solid harmonics R_{l,m} (homogeneous degree-l harmonic polynomials, nabla^2 R = 0) as
// Cartesian-monomial expansions, l = 0..3 (s,p,d,f), ordered m = -l..+l.  Only relative coefficients
// within each m matter (overall scale is normalised away by the evaluator); harmonicity + same-l
// orthogonality under the Gaussian inner product is what the unit test checks.
module;
#include <vector>
#include <cassert>
module qchem.BasisSet.Molecule.Evaluators.PG_Spherical_MnD.SolidHarmonics;

namespace qchem::BasisSet::Molecule::Evaluators::PG_Spherical_MnD
{
using P = Cart::Polarization;   // P(a,b,c) == monomial x^a y^b z^c

std::vector<std::vector<CartTerm>> SphericalShell(int l)
{
    switch (l)
    {
    case 0: return {
        { {P(0,0,0), 1.0} },                                              // s
    };
    case 1: return {
        { {P(0,1,0), 1.0} },                                             // m=-1  y
        { {P(0,0,1), 1.0} },                                             // m= 0  z
        { {P(1,0,0), 1.0} },                                             // m=+1  x
    };
    case 2: return {
        { {P(1,1,0), 1.0} },                                            // m=-2  xy
        { {P(0,1,1), 1.0} },                                            // m=-1  yz
        { {P(0,0,2), 1.0}, {P(2,0,0),-0.5}, {P(0,2,0),-0.5} },          // m= 0  z^2 - (x^2+y^2)/2
        { {P(1,0,1), 1.0} },                                            // m=+1  xz
        { {P(2,0,0), 1.0}, {P(0,2,0),-1.0} },                           // m=+2  x^2 - y^2
    };
    case 3: return {
        { {P(2,1,0), 3.0}, {P(0,3,0),-1.0} },                          // m=-3  3x^2 y - y^3
        { {P(1,1,1), 1.0} },                                           // m=-2  xyz
        { {P(0,1,2), 4.0}, {P(2,1,0),-1.0}, {P(0,3,0),-1.0} },         // m=-1  y(4z^2 - x^2 - y^2)
        { {P(0,0,3), 2.0}, {P(2,0,1),-3.0}, {P(0,2,1),-3.0} },         // m= 0  z(2z^2 - 3x^2 - 3y^2)
        { {P(1,0,2), 4.0}, {P(3,0,0),-1.0}, {P(1,2,0),-1.0} },         // m=+1  x(4z^2 - x^2 - y^2)
        { {P(2,0,1), 1.0}, {P(0,2,1),-1.0} },                          // m=+2  z(x^2 - y^2)
        { {P(3,0,0), 1.0}, {P(1,2,0),-3.0} },                          // m=+3  x(x^2 - 3y^2)
    };
    default:
        assert(false && "PG_Spherical_MnD::SolidHarmonics: only l=0..3 (s,p,d,f) supported");
        return {};
    }
}

} //namespace

// File: Common/CMath.C  The <cmath> leaf: std math functions + the cheap stable numeric primitives.
// (Imported by code INSIDE the math library; outside code imports the qchem.Math umbrella instead.)
module;
#include <cmath>
export module qchem.CMath;
export import qchem.IntPow;
export import qchem.Constants;
export import qchem.Factorials;

export
{
    using std::sqrt;
    using std::cbrt;
    using std::fabs;
    using std::abs;
    using std::floor;
    using std::ceil;
    using std::pow;
    using std::exp;
    using std::log;
    using std::log10;
    using std::sin;
    using std::cos;
    using std::acos;
    using std::erf;
    using std::erfc;
    using std::isfinite;
    using std::max;
    using std::min;
    using std::lround;

}


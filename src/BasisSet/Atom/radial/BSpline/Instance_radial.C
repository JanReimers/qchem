
module;
#include "radial/BSpline/GLQuadrature.H"

import Common.IntPow;
import oml;

module qchem.Basisset.Atom.radial.BSpline.Rk;

#include "./Imp/Rk.Ci"

namespace BSpline
{
    template class RkCache<6>;
    template class RkEngine<6>;
}
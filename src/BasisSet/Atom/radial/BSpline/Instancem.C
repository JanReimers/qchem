module;
#include "radial/BSpline/Rk.H"
#include "radial/BSpline/IBS_Common.H"
#include "radial/BSpline/IE_Primatives.H"
#include <iostream>
#include <cassert>
#include <bspline/operators/Derivative.h>
#include <iostream>


module qchem.BasisSet.Atom.l.BSplineBS;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.BasisSet.Imp.HeapDB;
import Common.IntPow;
import oml;

#include "../../l/Imp/BSpline_BS.Ci"
#include "../../l/Imp/BSpline_IBS.Ci"
#include "../../l/Imp/BSpline_BF.Ci"

namespace Atoml::BSpline
{
    //Initial attempts to loop this with meta programming failed.
    // template class BasisSet<3>;
    // template class BasisSet<4>;
    // template class BasisSet<5>;
    template class BasisSet<6>;
    template class Orbital_IBS<6>;
    template class BasisFunction<6>;

    // template class BasisSet<7>;
    // template class BasisSet<8>;
    // template class BasisSet<9>;
    // template class BasisSet<10>;
    // template class BasisSet<11>;
}


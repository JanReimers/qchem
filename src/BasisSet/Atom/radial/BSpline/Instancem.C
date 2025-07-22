module;
#include "radial/BSpline/IBS_Common.H"
#include "radial/BSpline/IE_Primatives.H"
#include <iostream>
#include <cassert>
#include <bspline/operators/Derivative.h>
#include <cmath>


module qchem.BasisSet.Atom.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.BasisSet.Imp.HeapDB;
import Common.IntPow;
import qchem.BasisSet.Imp.Cache4;
import oml;

#include "../../l/Imp/BSpline_BS.Ci"
#include "../../l/Imp/BSpline_IBS.Ci"
#include "../../l/Imp/BSpline_BF.Ci"
namespace Atoml::BSpline
{
    template class BasisSet<6>;
    template class Orbital_IBS<6>;
    template class BasisFunction<6>;
}


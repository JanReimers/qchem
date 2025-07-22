module;
#include <iostream>
#include <cassert>
#include <bspline/operators/Derivative.h>
#include <cmath>


module qchem.BasisSet.Atom.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.BS_Common; 
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.BasisSet.Imp.HeapDB;
import Common.IntPow;
import qchem.BasisSet.Imp.Cache4;
import qchem.Basisset.Atom.radial.BSpline.IE_Primatives;
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


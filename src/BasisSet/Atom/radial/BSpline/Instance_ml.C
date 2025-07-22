module;
#include <iostream>
#include <cassert>
#include <bspline/operators/Derivative.h>
#include <cmath>

#include "radial/BSpline/IBS_Common.H"
#include "radial/BSpline/BS_Common.H"

import qchem.Symmetry.Ylm;


module qchem.BasisSet.Atom.ml.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.BasisSet;
import qchem.Symmetry.Yl;
import qchem.BasisSet.Imp.HeapDB;
import Common.IntPow;
import qchem.Symmetry.AtomEC;
import qchem.Basisset.Atom.radial.BSpline.IE_Primatives;

import oml;

#include "../../ml/Imp/BSpline_BS.Ci"
#include "../../ml/Imp/BSpline_IBS.Ci"
#include "../../ml/Imp/BSpline_BF.Ci"

namespace Atom_ml::BSpline
{
    template class BasisSet<6>;
    template class Orbital_IBS<6>;
    template class BasisFunction<6>;
}

// File: Atom/radial/BSpline/Instance.C  Make instance of BSpline templates

#include <cmath>
import qchem.BasisSet.Imp.IEClient;

#include "BS_Common.Ci"
#include "IBS_Common.Ci"
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Atom.l.BSplineBS;
import qchem.BasisSet.Atom.ml.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Imp.IEClient;
import qchem.Basisset.Atom.radial.BSpline.IE_Primatives;

namespace BSpline
{
    template class IrrepBasisSet<6>;
    
    // template class RkEngine<6>;
    // template class IE_Fit<6>;
    // template class IE_DFT<double,6>;
    // template class IE_Overlap<double,6>;
    // template class IE_Inv_r1<double,6>;
    // template class IE_Kinetic<double,6>;
    template class BS_Common<6>;
    template class IE_Primatives<6>;
}

namespace Atoml::BSpline
{
    //Initial attempts to loop this with meta programming failed.
    // template class BasisSet<3>;
    // template class BasisSet<4>;
    // template class BasisSet<5>;
    template class BasisSet<6>;
    // template class BasisSet<7>;
    // template class BasisSet<8>;
    // template class BasisSet<9>;
    // template class BasisSet<10>;
    // template class BasisSet<11>;
}

namespace Atom_ml::BSpline
{
    //Initial attempts to loop this with meta programming failed.
    // template class BasisSet<3>;
    // template class BasisSet<4>;
    // template class BasisSet<5>;
    template class BasisSet<6>;
    // template class BasisSet<7>;
    // template class BasisSet<8>;
    // template class BasisSet<9>;
    // template class BasisSet<10>;
    // template class BasisSet<11>;
}

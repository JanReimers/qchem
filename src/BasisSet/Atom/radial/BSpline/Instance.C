// File: Atom/radial/BSpline/Instance.C  Make instance of BSpline templates

#include <cmath>
#include "radial/BSpline/IE_Primatives.H"
#include "BFGrouper.Ci"
#include "BS_Common.Ci"
#include "IBS_Common.Ci"
#include "IE_Primatives.Ci"
#include "IE.Ci"
#include "IEC.Ci"
#include "Rk.Ci"
#include "../../ml/BSpline_BS.Ci"
#include "../../ml/BSpline_IBS.Ci"
#include "../../ml/BSpline_BF.Ci"
import qchem.BasisSet.Atom.l.BSplineBS;
namespace BSpline
{
    template class IrrepBasisSet<6>;
    template class RkEngine<6>;
    template class IE_Fit<6>;
    template class BS_Common<6>;
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

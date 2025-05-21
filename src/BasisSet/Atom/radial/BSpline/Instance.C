// File: Atom/radial/BSpline/Instance.C  Make instance of BSpline templates

#include "BFGrouper.Ci"
#include "BS_Common.Ci"
#include "IBS_Common.Ci"
#include "IE_Primatives.Ci"
#include "IE.Ci"
#include "IEC.Ci"
#include "Rk.Ci"

namespace BSpline
{
template class BS_Common<6>;
template class BFGrouper<6>;
template class IrrepBasisSet<6>;
template class IE_Primatives<6>;
template class IE_Overlap<double,6>;
template class IE_Kinetic  <double,6>;
template class IE_Inv_r1<double,6>;
template class IE_DFT<double,6>;
template class IE_BS_2E<double,6>;
template class IE_Fit<6>;
template class IE_RKBL<double,6>;
template class IE_RKBS<double,6>;
}

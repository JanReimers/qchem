// File: BasisSet/Imp/Fit_Types.C  What a Vxc fit family forces (Forces; doc/Pins.md pin 30, Group face step 3b).
module;
#include <string>
#include <stdexcept>
module qchem.BasisSet.Fit_Types;
import qchem.Symmetry.Group;

namespace qchem::BasisSet
{

// THE ONE PLACE a Vxc fit kind's group is stated (controller 2026-11-21): every fit class's Forces() and the pin-30 pre-flight ask this
// switch.  If VxcFit ever grows a class per kind, this descriptor moves onto that class.
Outcome<Symmetry::GroupPtr,std::string> Forces(VxcFit fit, const Symmetry::Frame& frame)
{
    using R=Outcome<Symmetry::GroupPtr,std::string>;
    using namespace Symmetry;
    switch (fit)
    {
        case VxcFit::PlaneWave:
            if (!frame.lattice) return R::Fail("the plane-wave Vxc fit imposes H_R = T, and there is no lattice here for its ball -- To fix this, please state the Vxc "
                                               "fit Delta or Auxiliary");
            return R::Ok(TranslationGroup(*frame.lattice, kExactTolerance));
        case VxcFit::Auxiliary:
            if (frame.lattice) return R::Ok(TranslationGroup(*frame.lattice, kExactTolerance));   // lattice-summed at Gamma
            if (frame.centre)  return R::Ok(PointRotationGroup(*frame.centre, kExactTolerance));  // the l = 0 set about the site
            return R::Ok(TrivialGroup());                                                         // atom-centred Gaussians of every l
        case VxcFit::Delta:
        case VxcFit::Multipole:
        case VxcFit::Auto:
            return R::Ok(TrivialGroup());
    }
    throw std::logic_error("BasisSet::Forces: a Vxc fit family with no group");
}

} //namespace

// File: ChargeDensity.C  Interface for the charge density category.
module;
#include <memory>   // unique_ptr: densities are BUILT and handed over (V1.25)
#include <cassert>
module qchem.ChargeDensity;
import qchem.Symmetry.Spin;
import qchem.Blaze;

namespace qchem::ChargeDensity
{

// (The polarized CONTAINER -- tPolarized_CD, its Fourier/HF CRTP mixins and the Up/Down partner guard --
//  is gone: V1.37.  A polarized density is ONE tComposite_CD over full Irreps whose channels are VIEWS;
//  see Imp/CompositeCD.C, where the ↑/↓ summation tree those bodies expressed lives on as the
//  per-channel grouping.)

//----------------------------------------------------------------------------------
//
//  The magnetization m(r) = rho_up - rho_down.
//
template <class T> tSpinDensity<T>::tSpinDensity(std::unique_ptr<tDM_CD<T>> up,
                                                 std::unique_ptr<tDM_CD<T>> down)
: itsSpinUpCD  (std::move(up  ))
, itsSpinDownCD(std::move(down))
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};
// No destructor: the unique_ptr members free the channels, and being move-only they delete the copy
// operations whose absence made the raw-pointer form a double-delete (V1.25).

template <class T> double tSpinDensity<T>::operator()(const rvec3_t& r) const
{
    // No UT coverage
    return (*itsSpinUpCD)(r) - (*itsSpinDownCD)(r);
}

template <class T> rvec3_t tSpinDensity<T>::Gradient  (const rvec3_t& r) const
{
    // No UT coverage
    return itsSpinUpCD->Gradient(r) - itsSpinDownCD->Gradient(r);
}

template class tSpinDensity<double>;
template class tSpinDensity<dcmplx>;

} //namespace

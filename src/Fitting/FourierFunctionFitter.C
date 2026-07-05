// File: Fitting/FourierFunctionFitter.C  Plane-wave (Fourier) function fitter.
//
// The orthonormal/exact sibling of the molecular fitters, for the XC (overlap-metric) route: DoFit just
// RECEIVES the V-tilde the term computed (ForwardGrid on v_xc) and Overlap assembles <i|v_xc|j> directly
// (no kernel), delegating to the basis's reciprocal-space assembly (Band_FT_IBS).
//
// The DENSITY (Hartree) route moved to OrthoFunctionFitter (the core FunctionFitter_Density<dcmplx> face,
// reached through the factory MakeDensityFitter) -- PW_Hartree no longer stack-allocates this.  This XC-only
// remnant is retired when the XC/Scalar face is routed through its own ortho scalar fitter (the next
// increment); until then PW_XC still drives it concretely.
module;
#include <cassert>
#include <ostream>
export module qchem.Fitting.FourierFunctionFitter;
export import qchem.FourierMap;                // FourierMap (the G-space coefficients DoFit receives)
import qchem.Fitting.Types;                    // robs_t<T> (= Orbital_1E_IBS<T>, the common orbital base)
import qchem.BasisSet.Band_FT_IBS;             // the reciprocal-space contraction the fitter delegates to
import qchem.Blaze;                            // hmat_t<dcmplx>

export namespace qchem::Fitting
{

//! \brief The plane-wave XC "fitter": the potential's G-space coefficients (a FourierMap V-tilde) ARE the
//! fit (orthonormal, exact -- no metric solve), so DoFit stores them and Overlap delegates to Band_FT_IBS.
//! Standalone (no abstract base).  Driven concretely by PW_XC (the density/Hartree route is now
//! OrthoFunctionFitter behind the FunctionFitter_Density face).
class FourierFunctionFitter
{
public:
    //! The "fit": receive the potential's V-tilde (the term ran ForwardGrid on v_xc).  Orthonormal
    //! exactness => nothing to solve, just store; Overlap then reads the map with no kernel.
    void DoFit(const FourierMap& map) {itsMap=map;}

    //! XC matrix <i|v_xc|j>: assemble directly from the stored V-tilde (no kernel).
    hmat_t<dcmplx> Overlap(const robs_t<dcmplx>* bs) const
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "FourierFunctionFitter::Overlap requires a Band_FT_IBS (plane-wave) basis");
        return pw->Overlap(itsMap);
    }

    std::ostream& Write(std::ostream& os) const
        {return os << "FourierFunctionFitter (G-space map)" << std::endl;}

private:
    FourierMap itsMap;   //!< the fit = the potential's V-tilde (G-space coefficients).
};

} //namespace

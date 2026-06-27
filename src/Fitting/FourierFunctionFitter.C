// File: Fitting/FourierFunctionFitter.C  Plane-wave (Fourier) function fitter.
//
// The orthonormal/exact sibling of the molecular fitters.  For plane waves the density's rho-tilde IS the
// fit (the basis is orthonormal, so there is no metric solve): DoFit just RECEIVES the FourierMap the
// charge density already computed (MakeFourierDensity + BZ sum), and the contraction delegates to the
// basis's reciprocal-space assembly (Band_FT_IBS).  The orthonormal {G} basis makes the overlap/Coulomb
// metric distinction DEGENERATE -- the SAME map serves Hartree (Repulsion, with the 4pi/G^2 kernel) and XC
// (Overlap, no kernel) -- so this is ONE standalone class, not two abstract faces.
//
// It is used CONCRETELY (PW_Hartree/PW_XC stack-allocate it), never through a FunctionFitter_* base, so it
// carries NO abstract bases and NO NA-asserts: just the three real operations it actually performs.
module;
#include <cassert>
#include <ostream>
export module qchem.Fitting.FourierFunctionFitter;
export import qchem.FourierMap;                // FourierMap (the G-space coefficients DoFit receives)
import qchem.Fitting.Types;                    // obs_t<T> (= Orbital_1E_IBS<T>, the common orbital base)
import qchem.BasisSet.Band_FT_IBS;             // the reciprocal-space contraction the fitter delegates to
import qchem.Blaze;                            // hmat_t<dcmplx>

export namespace qchem::Fitting
{

//! \brief The plane-wave "fitter": the density's G-space coefficients (a FourierMap) ARE the fit
//! (orthonormal, exact -- no metric solve), so DoFit stores them and the contraction delegates to
//! Band_FT_IBS.  Standalone (no abstract base): the metric is degenerate on {G}, so one object serves both
//! Hartree (Repulsion) and XC (Overlap).  Driven concretely by PW_Hartree / PW_XC.
class FourierFunctionFitter
{
public:
    //! The "fit": receive the G-space coefficients -- the density's rho-tilde (Hartree) OR the potential's
    //! V-tilde (XC, the term ran ForwardGrid on v_xc).  Orthonormal exactness => nothing to solve, just store.
    //! Repulsion then reads the map with the 4pi/G^2 Coulomb kernel; Overlap reads it with no kernel.
    void DoFit(const FourierMap& map) {itsMap=map;}

    //! Coulomb (Hartree) matrix: delegate to the orbital basis's reciprocal-space Poisson assembly.  The
    //! orbital basis arrives as the common Orbital_1E_IBS base; cast down to the G-space capability.
    hmat_t<dcmplx> Repulsion(const obs_t<dcmplx>* bs) const
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "FourierFunctionFitter::Repulsion requires a Band_FT_IBS (plane-wave) basis");
        double Eh;
        return pw->Repulsion(itsMap, Eh);
    }

    //! XC matrix <i|v_xc|j>: assemble directly from the stored V-tilde (no kernel).
    hmat_t<dcmplx> Overlap(const obs_t<dcmplx>* bs) const
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "FourierFunctionFitter::Overlap requires a Band_FT_IBS (plane-wave) basis");
        return pw->Overlap(itsMap);
    }

    std::ostream& Write(std::ostream& os) const
        {return os << "FourierFunctionFitter (G-space map)" << std::endl;}

private:
    FourierMap itsMap;   //!< the fit = G-space coefficients: rho-tilde for Hartree, V-tilde for XC.
};

} //namespace

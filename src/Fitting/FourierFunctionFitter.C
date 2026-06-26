// File: Fitting/FourierFunctionFitter.C  Plane-wave (Fourier) function fitter.
//
// The orthonormal/exact sibling of FunctionFitterImp.  For plane waves the density's rho-tilde IS the fit
// (the basis is orthonormal, so there is no metric solve): DoFit just RECEIVES the FourierMap the charge
// density already computed (MakeFourierDensity + BZ sum), and the contraction delegates to the basis's
// reciprocal-space assembly (Band_FT_IBS::Repulsion).  Because it is a FunctionFitter<dcmplx>, a Hamiltonian
// term drives it exactly as FittedVee drives a FunctionFitter<double>:  DoFit(density) then Repulsion(obs).
//
// Only the Fourier (Hartree) path is wired today; the Gaussian client fits, the SCF-lifecycle helpers, and
// the scalar/XC route NA-assert -- the XC path is grid-based and does not (yet) flow through this fitter.
module;
#include <cassert>
#include <ostream>
export module qchem.Fitting.FourierFunctionFitter;
export import qchem.Fitting.FunctionFitter;   // FunctionFitter<dcmplx>, the *FFClients, FourierMap (re-exported)
import qchem.Fitting.Types;                    // obs_t<T> (= Orbital_1E_IBS<T>, the common orbital base)
import qchem.BasisSet.Band_FT_IBS;             // the reciprocal-space contraction the fitter delegates to
import qchem.Blaze;                            // hmat_t<dcmplx>, rvec3_t

export namespace qchem::Fitting
{

//! \brief The plane-wave "fitter": the density's G-space coefficients (a FourierMap) ARE the fit
//! (orthonormal, exact -- no metric solve), so DoFit stores them and the contraction delegates to
//! Band_FT_IBS.  A FunctionFitter<dcmplx>, so the Kohn-Sham terms drive it like the molecular fitter.
class FourierFunctionFitter : public virtual FunctionFitter<dcmplx>
{
public:
    //! The "fit": receive the G-space coefficients -- the density's rho-tilde (Hartree) OR the potential's
    //! V-tilde (XC, the term ran ForwardGrid on v_xc).  Orthonormal exactness => nothing to solve, just store.
    //! Repulsion then reads the map with the 4pi/G^2 Coulomb kernel; Overlap reads it with no kernel.
    virtual void DoFit(const ProjectedDensity_FT& map) override {itsMap=map;}

    //! Coulomb (Hartree) matrix: delegate to the orbital basis's reciprocal-space Poisson assembly.  The
    //! orbital basis arrives as the common Orbital_1E_IBS base; cast down to the G-space capability.
    virtual hmat_t<dcmplx> Repulsion(const obs_t<dcmplx>* bs) const override
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "FourierFunctionFitter::Repulsion requires a Band_FT_IBS (plane-wave) basis");
        double Eh;
        return pw->Repulsion(itsMap, Eh);
    }

    //! XC matrix <i|v_xc|j>: assemble directly from the stored V-tilde (no kernel).
    virtual hmat_t<dcmplx> Overlap(const obs_t<dcmplx>* bs) const override
    {
        auto pw=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
        assert(pw && "FourierFunctionFitter::Overlap requires a Band_FT_IBS (plane-wave) basis");
        return pw->Overlap(itsMap);
    }

    virtual std::ostream& Write(std::ostream& os) const override
        {return os << "FourierFunctionFitter (G-space map)" << std::endl;}

    // --- NA: the AO client fits + the SCF-lifecycle helpers don't flow through the plane-wave fitter. ---
    virtual void   DoFit(const ScalarFFClient&)  override
        {assert(false && "FourierFunctionFitter::DoFit(ScalarFFClient): PW v_xc arrives as a FourierMap (ForwardGrid)");}
    virtual void   DoFit(const ProjectedDensity_AO&) override
        {assert(false && "FourierFunctionFitter::DoFit(ProjectedDensity_AO): the PW density arrives as a FourierMap");}
    virtual void   ReScale(double) override
        {assert(false && "FourierFunctionFitter::ReScale: SCF mixing is done on the PW density, not the fitter");}
    virtual void   FitMixIn(const FunctionFitter<dcmplx>&,double) override
        {assert(false && "FourierFunctionFitter::FitMixIn: not used by the PW path");}
    virtual double FitGetChangeFrom(const FunctionFitter<dcmplx>&) const override
        {assert(false && "FourierFunctionFitter::FitGetChangeFrom: not used by the PW path"); return 0.0;}
    virtual double FitGetSelfRepulsion() const override
        {assert(false && "FourierFunctionFitter::FitGetSelfRepulsion: the term takes E_H via DM_Contract"); return 0.0;}
    virtual double Integral() const override
        {assert(false && "FourierFunctionFitter::Integral: not used by the PW path"); return 0.0;}
    virtual double  operator()(const rvec3_t&) const override
        {assert(false && "FourierFunctionFitter::operator(): not evaluated pointwise on the PW path"); return 0.0;}
    virtual rvec3_t Gradient  (const rvec3_t&) const override
        {assert(false && "FourierFunctionFitter::Gradient: not evaluated pointwise on the PW path"); return rvec3_t(0,0,0);}

private:
    FourierMap itsMap;   //!< the fit = G-space coefficients: rho-tilde for Hartree, V-tilde for XC.
};

} //namespace

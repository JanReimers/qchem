// File: Hamiltonian/Internal/GGAFunctional.C  A GGA exchange-correlation functional: its own input type, sigma REQUIRED.
module;
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
export module qchem.Hamiltonian.Internal.GGAFunctional;
import qchem.Streamable;
export import qchem.Types;   // rvec_t -- the batch records

export namespace qchem::Hamiltonian
{

//! \brief THE INPUT of a GGA evaluation, spin-native: \f$\rho_\uparrow,\rho_\downarrow\f$ and the \f$\sigma\f$ triple
//! \f$\sigma_{\sigma\tau}=\nabla\rho_\sigma\cdot\nabla\rho_\tau\f$, one entry per point in the caller's order.
//!
//! THE GRADIENT IS A TYPE, NOT A RUN-TIME PROPERTY (user, 2026-11-20: "if GGA/PBE needs a gradient of object X, grad op is either
//! 1) In the X interface or b) not in the X interface.  If b the compile fails").  \f$\sigma\f$ is not optional here and there is no
//! "has gradient" query: the only constructor takes all five arrays, so code that cannot supply \f$\nabla\rho\f$ cannot make one, and
//! an LDA's input is a different type (\c LDAInput).  An unpolarized run is the term's \f$\rho/2\f$, \f$\sigma/4\f$ collapse (pin 5).
struct GGAInput
{
    //! THROWS if the five arrays are not one length (a caller's bookkeeping slip, caught at the record, not inside libxc).
    GGAInput(rvec_t up_, rvec_t dn_, rvec_t sigUU_, rvec_t sigUD_, rvec_t sigDD_)
        : up(std::move(up_)), dn(std::move(dn_)), sigUU(std::move(sigUU_)), sigUD(std::move(sigUD_)), sigDD(std::move(sigDD_))
    {
        const size_t n=up.size();
        if (dn.size()!=n || sigUU.size()!=n || sigUD.size()!=n || sigDD.size()!=n)
            throw std::invalid_argument("GGAInput: rho_up, rho_dn and the sigma triple must have one entry per point ("+std::to_string(n)+")");
    }
    size_t size() const {return up.size();}

    rvec_t up, dn;                 //!< \f$\rho_\uparrow,\rho_\downarrow\f$
    rvec_t sigUU, sigUD, sigDD;    //!< \f$\sigma_{\uparrow\uparrow},\sigma_{\uparrow\downarrow},\sigma_{\downarrow\downarrow}\f$
};

//! \brief THE OUTPUT of a GGA: \f$e\f$ per unit volume, \f$v_{\rho\sigma}=\partial e/\partial\rho_\sigma\f$ and
//! \f$v_{\sigma\sigma\tau}=\partial e/\partial\sigma_{\sigma\tau}\f$ (libxc's \c vsigma layout).  A composite SUMS outputs.
struct GGAOutput
{
    //! Zeroed, one entry per point of \a in: what an \c Evaluate accumulates into.
    explicit GGAOutput(const GGAInput& in)
        : e(in.size(),0.0), vUp(in.size(),0.0), vDn(in.size(),0.0), vsUU(in.size(),0.0), vsUD(in.size(),0.0), vsDD(in.size(),0.0) {}

    rvec_t e;                      //!< energy density per unit volume
    rvec_t vUp, vDn;               //!< \f$\partial e/\partial\rho_\uparrow,\ \partial e/\partial\rho_\downarrow\f$
    rvec_t vsUU, vsUD, vsDD;       //!< \f$\partial e/\partial\sigma_{\uparrow\uparrow},\ \partial e/\partial\sigma_{\uparrow\downarrow},\ \partial e/\partial\sigma_{\downarrow\downarrow}\f$
};

//! \brief A GGA exchange-correlation functional: ONE face, the batch over a \c GGAInput (GGA design note §2, revised 2026-11-20).
//!
//! It is deliberately NOT an \c ExFunctional: that face's scalar \f$\rho\mapsto v\f$ cannot be answered without \f$\sigma\f$, and a
//! face a class can only throw from is a run-time trap.  So no \f$\rho\f$-only call compiles on a GGA, and an XC term that takes
//! a \c GGAFunctional must also be handed something that produces \f$\nabla\rho\f$ (the term's constructor, stage 2).  Data-free.
class GGAFunctional
    : public virtual Streamable
{
public:
    //! ADDS this functional's \f$e\f$, \f$v_\rho\f$, \f$v_\sigma\f$ to \a out (one entry per point of \a in), so a composite is a loop.
    virtual void Evaluate(const GGAInput& in, GGAOutput& out) const=0;
};

//! \brief A SUM OF GGA PARTS behaving as one (PBE = GGA_X_PBE + GGA_C_PBE): the GGA twin of \c CompositeExFunctional, for the same
//! reason (one gather of the summed potential, not one per part).
class CompositeGGAFunctional
    : public virtual GGAFunctional
{
public:
    //! \a parts must be non-empty (THROWS otherwise: a sum of no functionals is not a functional).
    explicit CompositeGGAFunctional(std::vector<std::shared_ptr<const GGAFunctional>> parts) : itsParts(std::move(parts))
    {
        if (itsParts.empty()) throw std::invalid_argument("CompositeGGAFunctional: a sum of no functionals is not a functional");
    }
    virtual void Evaluate(const GGAInput& in, GGAOutput& out) const {for (const auto& p : itsParts) p->Evaluate(in,out);}
    virtual std::ostream& Write(std::ostream& os) const {for (const auto& p : itsParts) p->Write(os); return os;}

private:
    std::vector<std::shared_ptr<const GGAFunctional>> itsParts;
};

} //namespace

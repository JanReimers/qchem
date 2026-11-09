// File: BasisSet/CoulombDeltaFit_IBS.C  A δ basis WITH its Coulomb integrals: the point set a fitted Hartree term runs on.
//
// V-ONE-OF-N 3b-ii (user ruling 2026-11-08 (f), doc/Records/VOneOfN_D2Design.md §1a): 1/r12 plays two roles, and on a mesh the
// PHYSICS one is the point set's Coulomb integrals.  A Gaussian fit basis answers it analytically (⟨ab|1/r|f⟩, ⟨f|1/r|g⟩); a plane-
// wave basis as k(G)⟨ab|G⟩; a δ basis through a CoulombOperator (direct sum, radial Poisson, multipole).  So the operator belongs
// to the BASIS, not to the fitter that runs on it, and the fitter asks the basis for its potential.
//
// A DECORATOR OF THE δ BASIS, as MultipoleFit_IBS is: the same functions, the same forward and adjoint, plus the integrals.  The δ
// basis itself stays operator-free because the XC quadrature uses the same point set and needs no 1/r12.
module;
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
export module qchem.BasisSet.CoulombDeltaFit_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.Mesh.CoulombOperator;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Symmetry;
import qchem.Types;
export namespace qchem::BasisSet
{

//! \brief A point set's COULOMB INTEGRALS, asked of the basis (the physics round of a fitted Hartree term).  \a q are the projections
//! \f$\langle f_a|\rho\rangle\f$ on the basis's functions; the answers are on the same functions.  Data-free.
class PointCoulombIntegrals
{
public:
    virtual ~PointCoulombIntegrals() = default;
    virtual rvec_t Potential    (const rvec_t& q) const=0;                               //!< the potential's coefficients
    virtual double SelfRepulsion(const rvec_t& q, const rvec_t& v) const=0;              //!< \f$\langle\tilde\rho|\tilde\rho\rangle\f$
    virtual double HartreeEnergy(double trDV, const rvec_t& q, const rvec_t& v) const=0; //!< the first-order \f$E_H\f$
    virtual bool   UsesTrDV     () const=0;                                              //!< does \c HartreeEnergy read \a trDV?
};

//! \brief The δ basis \a delta with the Coulomb operator \a op over its points.
class CoulombDeltaFit_IBS
    : public virtual rFIT_SF_ABS
    , public virtual PointCoulombIntegrals
    , public         IrrepBasisSetImp<double>
    , public         Symmetry::IdentityRepresentation<double>
{
public:
    CoulombDeltaFit_IBS(std::shared_ptr<const rFIT_SF_ABS> delta, std::unique_ptr<qcMesh::CoulombOperator> op, const sym_t& sym)
        : IrrepBasisSetImp<double>(sym), itsDelta(std::move(delta)), itsOp(std::move(op))
    {
        if (!itsDelta || !itsOp) throw std::logic_error("CoulombDeltaFit_IBS: needs a δ basis and its Coulomb operator");
    }
    // --- the δ basis, unchanged
    bool   isOrtho        () const override {return itsDelta->isOrtho();}
    rvec_t Overlap        (const ScalarFunction<double>& f) const override {return itsDelta->Overlap(f);}
    rvec_t OverlapDiagonal() const override {return itsDelta->OverlapDiagonal();}
    rvec_t Charge         () const override {return itsDelta->Charge();}
    const PairProjection<double>& Overlap3C(const Orbital_DFT_IBS<double,double>& orb) const override {return itsDelta->Overlap3C(orb);}
    size_t        GetNumFunctions() const override {return itsDelta->GetNumFunctions();}
    std::string   Name           () const override {return itsDelta->Name();}
    //! The SAME basis as the δ one (the integrals do not change its functions), so the same ID and the same cached tables.
    std::string   BasisSetID     () const override {return itsDelta->BasisSetID();}
    std::ostream& Write          (std::ostream& os) const override {return itsDelta->Write(os);}
    // --- its Coulomb integrals
    rvec_t Potential    (const rvec_t& q) const override                               {return itsOp->Potential(q);}
    double SelfRepulsion(const rvec_t& q, const rvec_t& v) const override              {return itsOp->SelfRepulsion(q,v);}
    double HartreeEnergy(double trDV, const rvec_t& q, const rvec_t& v) const override {return itsOp->HartreeEnergy(trDV,q,v);}
    bool   UsesTrDV     () const override                                              {return itsOp->UsesTrDV();}

private:
    std::shared_ptr<const rFIT_SF_ABS>       itsDelta;
    std::unique_ptr<qcMesh::CoulombOperator> itsOp;
};

} //namespace

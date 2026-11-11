// File: BasisSet/Gaussian/Point/IrrepBasisSet.C  Molecule-generic, evaluator-templated orbital IBS mixins.
//
// The molecular analog of BasisSet/Radial/IrrepBasisSet.C (module qchem.BasisSet.Radial.IBS): the i,j /
// 3-centre / 4-centre matrix-build loops are basis-agnostic and driven purely by the evaluator's inline
// kernels, so they live here once, templated on the evaluator E (constrained by the Evaluators concepts),
// and every concrete molecular basis set (PolarizedGaussian today, SphericalGaussian / PRISM / libcint
// later) reuses them by instantiating the mixins with its own evaluator.
//
// As on the atom side, each mixin `dynamic_cast<const E&>(*this)`s to reach the evaluator base subobject
// of the final IBS (the IBS IS-A E -- see the view->base-subobject collapse in the PG tree).  The loop
// bodies are exactly the ones that used to be hand-written in PG_Cart/Imp/Orbital_IBS.C.
//
// Divergence from atom (expected -- "copy, don't share" until the common subset is proven, plan Goal C):
//  - 1E Kinetic here is the full Cartesian <p^2>=<-nabla^2> (no centrifugal split; that is atom-only).
//  - HF Direct/Exchange use the multi-centre M&D FourC kernel, NOT the atom Cache4/Grouper/Ak design.
module;
#include <cmath>
#include <concepts>
#include <stdexcept>
#include <string>
#include <cassert>
#include <memory>
#include <mutex>
#include <utility>
#include <vector>
export module qchem.BasisSet.Gaussian.Point.IBS;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Orbital_PP_IBS;            // the species-field service the PP mixin answers
import qchem.Math;                                      // Pi
import qchem.Math.SolidHarmonics;                       // CartGaussian, MultiplyR2, YlmCartesian
import qchem.Math.Angular;                               // CartTerm, Monomial
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Gaussian.Evaluators;      // concepts + generic 1E matrix builders
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD;   // NR_Evaluator: the evaluator whose blocks give libcint's derivative integrals
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;    // the derivative integrals (int1e_ip*), gated in tests/M_LibCint.C
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;   // PGData: the layout libcint's evaluator is built on
export import qchem.BasisSet.AoShellSource;              // the GetAoShells face (+ Symmetry::Molecule::AoShell)
export import qchem.BasisSet.BareCoulombSource;          // the BareCoulomb face (+ ERI4Block) -- ACBN0 (+U increment 3)
import qchem.Structure;
import qchem.Symmetry.Representation;               // DenseRepresentation (AoShellRepresentation, S3)
import qchem.Symmetry.Molecule.OperationRep;         // BuildOperationRep (the Cartesian SymOp overload)
import qchem.Types;

export namespace qchem::BasisSet::Gaussian
{

// --- The molecular orbital 1E IBS interface -------------------------------------------------------
// The non-templated interface that ADDS molecule-specific operations to the generic Real_OIBS
// (::qchem::BasisSet::Orbital_1E_IBS<double>) -- currently just GetAoShells (the point-group SALC seam).
// A client (PG::SymmetryAdapt) can Iterate<Gaussian::Orbital_1E_IBS> to reach every molecular orbital basis
// polymorphically, WITHOUT knowing its evaluator; the evaluator-templated EOrbital_1E_IBS<E> below supplies
// the evaluator-driven integral implementations.  Real-space by nature, so it lives on the molecule side and
// is absent from the dcmplx plane-wave path.
class Orbital_1E_IBS
    : public virtual ::qchem::BasisSet::Orbital_1E_IBS<double>
    , public virtual ::qchem::BasisSet::AoShellSource   // GetAoShells -- the face moved UP (2026-09-19) so a
                                                        // Bloch block over this one can answer it too
{
};

//! \brief A molecular Gaussian AO block's REPRESENTATION (V-BS-SYMMETRY S3): \f$M_g\f$ is the shell permutation times each
//! shell's angular rep, built from the block's own shells (\c GetAoShells) for an op \f$\{W|\tau\}\f$ in the Cartesian frame
//! -- the molecular twin of a GPW block at Γ.  Dense (n_AO is small); built per op on first use and cached by the base.
//! THROWS for an op that maps a shell onto no shell (not a symmetry of the molecule's basis).
class AoShellRepresentation
    : public virtual Symmetry::DenseRepresentation<double>
    , public virtual ::qchem::BasisSet::AoShellSource
{
protected:
    virtual rmat_t MakeRepMatrix(const Symmetry::SymOp& g) const override
    {
        return Symmetry::Molecule::BuildOperationRep(GetAoShells(), g, 1e-6);
    }
};

// --- 1E: Overlap / Kinetic(<p^2>) / Nuclear -------------------------------------------------------
// Two evaluator granularities, one mixin: if E delivers whole matrices (isM_1E_Evaluator -- an opaque
// assembler / block library), FORWARD; otherwise (is1E_Evaluator -- cheap per-element kernels) run the i,j
// loop.  The concrete IBS just instantiates with its evaluator; the framework adapts to whichever concept
// it satisfies.  (An evaluator that satisfied both would take the forward path, but the two are disjoint in
// practice -- a kernel evaluator wants the loop + shared cache, an assembler hands us the matrix.)
template <class E> requires (Evaluators::is1E_Evaluator<E> || Evaluators::isM_1E_Evaluator<E>)
class EOrbital_1E_IBS
    : public virtual Orbital_1E_IBS         // the molecule 1E interface (IS-A Real_OIBS + GetAoShells)
    , public virtual AoShellRepresentation  // every molecular Gaussian AO block states its representation (S3)
{
protected:
    virtual rsmat_t MakeOverlap() const
    {
        if constexpr (Evaluators::isM_1E_Evaluator<E>) return Cast().OverlapMatrix();
        else { const E& e=Cast(); rsmat_t S(e.size());
               for (auto i:e.indices()) for (auto j:e.indices(i)) S(i,j)=e.Overlap(i,j); return S; }
    }
    // <p^2>=<-nabla^2> building block: NO 1/2 (the Hamiltonian applies it) and NO centrifugal term (the
    // molecular Grad2 is already the full Cartesian -nabla^2).  See BasisSet/Orbital_1E_IBS.C.
    virtual rsmat_t MakeKinetic() const
    {
        if constexpr (Evaluators::isM_1E_Evaluator<E>) return Cast().KineticMatrix();
        else { const E& e=Cast(); rsmat_t S(e.size());
               for (auto i:e.indices()) for (auto j:e.indices(i)) S(i,j)=e.Grad2(i,j); return S; }
    }
    virtual rsmat_t MakeNuclear(const Structure* cl) const
    {
        if constexpr (Evaluators::isM_1E_Evaluator<E>) return Cast().NuclearMatrix(cl);
        else { const E& e=Cast(); rsmat_t S(e.size());
               for (auto i:e.indices()) for (auto j:e.indices(i)) S(i,j)=e.Nuclear(i,j,cl); return S; }
    }
    const E& Cast() const {return dynamic_cast<const E&>(*this);}

public:
    // d/dR (GeometryOptimizationPlan §13.2): libcint's derivative integrals on this block's own PGData, for the Cartesian PG_Cart evaluator
    // (increment 1).  Every other Gaussian evaluator refuses until it has them.
    virtual bool CanPullbackOverlap() const override {return kLibcint;}
    virtual bool CanPullbackKinetic() const override {return kLibcint;}
    virtual bool CanPullbackNuclear() const override {return kLibcint;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<double>& W, const Structure* st) const override
    {
        if constexpr (kLibcint) return Outcome<rvec3vec_t,AD::NoGradient>::Ok(Contract(LibCint<E>(st).OverlapGradient(), W));
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", "overlap pullback");
    }
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeKineticPullback(const hmat_t<double>& D, const Structure* st) const override
    {
        if constexpr (kLibcint) return Outcome<rvec3vec_t,AD::NoGradient>::Ok(Contract(LibCint<E>(st).KineticGradient(), D));
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", "kinetic pullback");
    }
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeNuclearPullback(const hmat_t<double>& D, const Structure* st) const override
    {
        if constexpr (kLibcint) return Outcome<rvec3vec_t,AD::NoGradient>::Ok(Contract(LibCint<E>(st).NuclearGradient(), D));
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", "nuclear pullback");
    }
private:
    static constexpr bool kLibcint=std::same_as<E, Evaluators::PG_Cart_MnD::NR_Evaluator>;
    //! A TEMPLATE, so it is instantiated only inside the kLibcint branch: for any other evaluator Cast() is no PGData.
    template <class EE> Evaluators::PG_LibCint::NR_Evaluator LibCint(const Structure* st) const
    {
        if (!st) throw std::invalid_argument("EOrbital_1E_IBS: a geometry pullback needs the structure");
        const EE& e=dynamic_cast<const EE&>(*this);
        return Evaluators::PG_LibCint::NR_Evaluator(static_cast<const Evaluators::PG_Cart_MnD::PGData&>(e), st);
    }
    //! Σ_ij A_ij dM_ij/dR_atom for every atom: a per-atom gradient-matrix set contracted with the adjoint.
    static rvec3vec_t Contract(const Evaluators::PG_LibCint::NR_Evaluator::GradientMatrices& dM, const hmat_t<double>& A)
    {
        rvec3vec_t F(dM.size(), rvec3_t(0,0,0));
        for (size_t at=0;at<dM.size();at++)
        {
            double x[3]={0,0,0};
            for (int k=0;k<3;k++)
            {
                const rmat_t& M=dM[at][k];
                if (M.rows()!=A.rows() || M.columns()!=A.columns())
                    throw std::invalid_argument("EOrbital_1E_IBS: the adjoint is "+std::to_string(A.rows())+" x "+std::to_string(A.columns())+
                                                ", the block's matrices are "+std::to_string(M.rows())+" x "+std::to_string(M.columns()));
                for (size_t i=0;i<M.rows();i++) for (size_t j=0;j<M.columns();j++) x[k]+=A(i,j)*M(i,j);
            }
            F[at]=rvec3_t(x[0],x[1],x[2]);
        }
        return F;
    }
};

// --- Species fields: the analytic local PP and KB projectors -----------------------------------------
// The kernels an evaluator supplies for a closed-form field.  OverlapWith is <i|g>, OperatorThreeC is
// <i|g|j>, and RepulsionThreeC is (ij|e^{-alpha|r-c|^2}), the Coulomb potential of a Gaussian charge.
template <class E> concept isPP_Evaluator = std::derived_from<E, Evaluators::Evaluator>
    && requires (const E e, size_t i, size_t j, const Math::CartGaussian& g, const rvec3_t& c, double a)
{
    {e.OverlapWith    (i, g)      } -> std::same_as<double>;
    {e.OperatorThreeC (i, j, g)   } -> std::same_as<double>;
    {e.RepulsionThreeC(i, j, c, a)} -> std::same_as<double>;
};

//! \brief The finite-structure \c Orbital_PP_IBS of a Gaussian basis, in closed form.  The same forms GPW
//! assembles under its lattice sums: a Gaussian short part, an erf long part, Gaussian-\f$Y_{lm}\f$
//! projectors.  Exact, so no mesh enters these integrals.  A field without the closed-form faces cannot be
//! assembled here, and \c CanAssembleField says so.
template <class E> requires isPP_Evaluator<E>
class EOrbital_PP_IBS : public virtual Orbital_PP_IBS<double>
{
public:
    virtual bool CanAssembleField(const SpeciesRadialField& f) const override
    {
        return dynamic_cast<const SpeciesRadialField_Gaussian*>(&f)
            && dynamic_cast<const SpeciesRadialField_GaussianCharge*>(&f);
    }
    virtual bool CanAssembleProjectors(const SpeciesProjectorSet& ps) const override
    {
        return dynamic_cast<const SpeciesProjectorSet_Gaussian*>(&ps);
    }
    virtual rsmat_t MakeSpeciesFieldMatrix(const Structure* st, const SpeciesRadialField& f, FieldRange range) const override
    {
        const E& e=dynamic_cast<const E&>(*this);
        const auto& gs=dynamic_cast<const SpeciesRadialField_Gaussian&>(f);
        const auto& gc=dynamic_cast<const SpeciesRadialField_GaussianCharge&>(f);
        rsmat_t V(e.size());
        for (const Site& a : st->Sites())
        {
            if (range!=FieldRange::Long)
                for (const auto& t : gs.AsGaussians(a.Z,FieldRange::Short))
                {
                    Math::CartGaussian g{a.R, t.alpha, Math::MultiplyR2({Math::CartTerm{Math::Monomial{0,0,0},t.c}}, t.n)};
                    for (auto i:e.indices()) for (auto j:e.indices(i)) V(i,j)+=e.OperatorThreeC(i,j,g);
                }
            if (range!=FieldRange::Short)
            {
                const GaussianCharge q=gc.AsGaussianCharge(a.Z);
                const double scale=-q.q*std::pow(q.alpha/Pi,1.5);
                for (auto i:e.indices()) for (auto j:e.indices(i)) V(i,j)+=scale*e.RepulsionThreeC(i,j,a.R,q.alpha);
            }
        }
        return V;
    }
    virtual rsmat_t MakeProjectorMatrix(const Structure* st, const SpeciesProjectorSet& ps) const override
    {
        const E& e=dynamic_cast<const E&>(*this);
        const auto& gp=dynamic_cast<const SpeciesProjectorSet_Gaussian&>(ps);
        const size_t n=e.size();
        rsmat_t V(n);
        for (const Site& a : st->Sites())
            for (size_t p=0; p<ps.Count(a.Z); p++)
            {
                const int l=ps.L(a.Z,p);
                const double D=ps.Weight(a.Z,p);
                const auto radial=gp.AsGaussians(a.Z,p);
                for (int m=-l; m<=l; m++)
                {
                    rvec_t b(n,0.0);
                    for (const auto& t : radial)                // c r^{2n} e^{-alpha r^2} x [r^l Y_lm]
                    {
                        Math::CartGaussian g{a.R, t.alpha, Math::MultiplyR2(Math::YlmCartesian(l,m), t.n)};
                        for (auto& c : g.terms) c.c*=t.c;
                        for (auto i:e.indices()) b[i]+=e.OverlapWith(i,g);
                    }
                    for (auto i:e.indices()) for (auto j:e.indices(i)) V(i,j)+=D*b[i]*b[j];
                }
            }
        return V;
    }

    // d/dR (GeometryOptimizationPlan §13.2): libcint's derivative integrals on this block's own PGData, each field Gaussian an auxiliary
    // shell (PG_LibCint::NR_Evaluator's PP gradients), for the Cartesian PG_Cart evaluator; every other evaluator refuses.  The VALUES
    // the projector pullback weights by (b_i = <i|β Y_lm>) are this engine's, exactly as MakeProjectorMatrix forms them.
    virtual bool CanPullbackSpeciesField(const SpeciesRadialField& f) const override {return kLibcint && CanAssembleField(f);}
    virtual bool CanPullbackProjectors(const SpeciesProjectorSet& ps) const override {return kLibcint && CanAssembleProjectors(ps);}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeSpeciesFieldPullback(const Structure* st, const SpeciesRadialField& f, FieldRange range,
                                                                        const hmat_t<double>& D) const override
    {
        if constexpr (kLibcint)
        {
            if (!CanAssembleField(f)) return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block, a field with no closed form", "species-field pullback");
            const auto lc=LibCintOf(st);
            const rmat_t Dd=Dense(D, lc.size());
            const auto& gs=dynamic_cast<const SpeciesRadialField_Gaussian&>(f);
            const auto& gc=dynamic_cast<const SpeciesRadialField_GaussianCharge&>(f);
            rvec3vec_t F(st->GetNumAtoms(), rvec3_t(0,0,0));
            auto add=[&](const std::vector<rvec3_t>& g, double s){ for (size_t A=0;A<g.size();A++) F[A]+=s*g[A]; };
            for (size_t A=0;A<st->GetNumAtoms();A++)
            {
                const Site& a=st->Sites()[A];
                if (range!=FieldRange::Long)
                    for (const auto& t : gs.AsGaussians(a.Z,FieldRange::Short))
                        add(lc.OperatorThreeCGradient(Math::CartGaussian{a.R, t.alpha, Math::MultiplyR2({Math::CartTerm{Math::Monomial{0,0,0},t.c}}, t.n)}, A, Dd), 1.0);
                if (range!=FieldRange::Short)
                {
                    const GaussianCharge q=gc.AsGaussianCharge(a.Z);
                    add(lc.GaussianRepulsionGradient(A, q.alpha, Dd), -q.q*std::pow(q.alpha/Pi,1.5));
                }
            }
            return Outcome<rvec3vec_t,AD::NoGradient>::Ok(F);
        }
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", "species-field pullback");
    }
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeProjectorPullback(const Structure* st, const SpeciesProjectorSet& ps, const hmat_t<double>& D) const override
    {
        if constexpr (kLibcint)
        {
            if (!CanAssembleProjectors(ps)) return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block, projectors with no closed form", "projector pullback");
            const E& e=dynamic_cast<const E&>(*this);
            const auto lc=LibCintOf(st);
            const size_t n=e.size();
            if (D.rows()!=n) throw std::invalid_argument("Gaussian EOrbital_PP_IBS::MakeProjectorPullback: D is "+std::to_string(D.rows())+" square for "+std::to_string(n)+" functions");
            const auto& gp=dynamic_cast<const SpeciesProjectorSet_Gaussian&>(ps);
            rvec3vec_t F(st->GetNumAtoms(), rvec3_t(0,0,0));
            // E = Σ_{a,p,m} w_p bᵀ D b, b_i = Σ_t <i|g_t>: so F[A] = Σ_{a,p,m} 2 w_p Σ_t Σ_i (D b)_i ∂<i|g_t>/∂R_A.
            for (size_t A=0;A<st->GetNumAtoms();A++)
            {
                const Site& a=st->Sites()[A];
                for (size_t p=0; p<ps.Count(a.Z); p++)
                {
                    const int l=ps.L(a.Z,p);
                    const double wp=ps.Weight(a.Z,p);
                    const auto radial=gp.AsGaussians(a.Z,p);
                    for (int m=-l; m<=l; m++)
                    {
                        std::vector<Math::CartGaussian> gs;
                        rvec_t b(n,0.0);
                        for (const auto& t : radial)
                        {
                            Math::CartGaussian g{a.R, t.alpha, Math::MultiplyR2(Math::YlmCartesian(l,m), t.n)};
                            for (auto& c : g.terms) c.c*=t.c;
                            for (auto i:e.indices()) b[i]+=e.OverlapWith(i,g);
                            gs.push_back(std::move(g));
                        }
                        rvec_t w(n,0.0);
                        for (size_t i=0;i<n;i++) { double s=0.0; for (size_t j=0;j<n;j++) s+=D(i,j)*b[j]; w[i]=2.0*wp*s; }
                        for (const auto& g : gs)
                        {
                            const auto dg=lc.OverlapWithGradient(g, A, w);
                            for (size_t B=0;B<dg.size();B++) F[B]+=dg[B];
                        }
                    }
                }
            }
            return Outcome<rvec3vec_t,AD::NoGradient>::Ok(F);
        }
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", "projector pullback");
    }
private:
    static constexpr bool kLibcint=std::same_as<E, Evaluators::PG_Cart_MnD::NR_Evaluator>;
    //! The libcint evaluator on this block's own PGData (a TEMPLATE body, instantiated only in the kLibcint branch).
    template <class EE=E> Evaluators::PG_LibCint::NR_Evaluator LibCintOf(const Structure* st) const
    {
        if (!st) throw std::invalid_argument("Gaussian EOrbital_PP_IBS: a geometry pullback needs the structure");
        return Evaluators::PG_LibCint::NR_Evaluator(static_cast<const Evaluators::PG_Cart_MnD::PGData&>(dynamic_cast<const EE&>(*this)), st);
    }
    static rmat_t Dense(const hmat_t<double>& D, size_t n)
    {
        if (D.rows()!=n) throw std::invalid_argument("Gaussian EOrbital_PP_IBS: D is "+std::to_string(D.rows())+" square for "+std::to_string(n)+" functions");
        rmat_t Dd(n,n); for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) Dd(i,j)=D(i,j);
        return Dd;
    }
};

// --- 3-centre (DFT fit): <ab|c> for each fit component c ------------------------------------------
// Forward Cast().*ThreeC_Matrix(fit) when E delivers the dense PairProjection (isM_DFT_Evaluator); else run the Make3C loop.
template <class E> requires (Evaluators::isDFT_Evaluator<E> || Evaluators::isM_DFT_Evaluator<E>)
class Orbital_DFT_IBS
    : public virtual ::qchem::BasisSet::Orbital_DFT_IBS<double>
{
protected:
    virtual PairProjection<double> MakeOverlap3C  (const rFIT_SF_ABS& c) const
    {
        if constexpr (Evaluators::isM_DFT_Evaluator<E>)
            return dynamic_cast<const E&>(*this).OverlapThreeC_Matrix(dynamic_cast<const E&>(c));
        else return Make3C(c, [](const E& aE, size_t ia, size_t ib, const E& cE, size_t ic)
                                   {return aE.OverlapThreeC(ia, aE, ib, cE, ic);});
    }
    virtual PairProjection<double> MakeRepulsion3C(const rFIT_CD_ABS& c) const
    {
        if constexpr (Evaluators::isM_DFT_Evaluator<E>)
            return dynamic_cast<const E&>(*this).RepulsionThreeC_Matrix(dynamic_cast<const E&>(c));
        else return Make3C(c, [](const E& aE, size_t ia, size_t ib, const E& cE, size_t ic)
                                   {return aE.RepulsionThreeC(ia, aE, ib, cE, ic);});
    }
public:
    // d/dR (GeometryOptimizationPlan §13.2): libcint's int3c2e_ip1/_ip2 on the two blocks' own PGData, for the Cartesian PG_Cart evaluator
    // when the fit block is one too (the same E, as MakeRepulsion3C casts it).  Every other pairing refuses until it has them.
    virtual bool CanPullbackRepulsion3C(const rFIT_CD_ABS& c) const override
    {
        if constexpr (kLibcint) return dynamic_cast<const E*>(&c)!=nullptr;
        else return false;
    }
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsion3CPullback(const rFIT_CD_ABS& c, const hmat_t<double>& D, const rvec_t& a,
                                                                       const Structure* st) const override
    {
        if constexpr (kLibcint)
        {
            if (!CanPullbackRepulsion3C(c)) return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block against a fit of another family", "3-centre repulsion pullback");
            if (!st) throw std::invalid_argument("Gaussian Orbital_DFT_IBS: a geometry pullback needs the structure");
            using LC=Evaluators::PG_LibCint::NR_Evaluator;
            using PGD=Evaluators::PG_Cart_MnD::PGData;
            const LC orb(static_cast<const PGD&>(dynamic_cast<const E&>(*this)), st), fit(static_cast<const PGD&>(dynamic_cast<const E&>(c)), st);
            if (D.rows()!=orb.size() || a.size()!=fit.size())
                throw std::invalid_argument("Gaussian Orbital_DFT_IBS::MakeRepulsion3CPullback: D is "+std::to_string(D.rows())+" square and a has "
                                            +std::to_string(a.size())+", for "+std::to_string(orb.size())+" orbital and "+std::to_string(fit.size())+" fit functions");
            rmat_t Dd(D.rows(), D.columns());
            for (size_t i=0;i<D.rows();i++) for (size_t j=0;j<D.columns();j++) Dd(i,j)=D(i,j);
            const std::vector<rvec3_t> g=orb.RepulsionThreeC_Gradient(fit, Dd, a);
            rvec3vec_t F(g.size()); for (size_t A=0;A<g.size();A++) F[A]=g[A];
            return Outcome<rvec3vec_t,AD::NoGradient>::Ok(F);
        }
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", "3-centre repulsion pullback");
    }
private:
    static constexpr bool kLibcint=std::same_as<E, Evaluators::PG_Cart_MnD::NR_Evaluator>;
    // For each fit component ic, build the symmetric (ia,ib) block via the supplied named 3-centre kernel
    // (which folds in all three normalizations).
    template <class Kernel>
    PairProjection<double> Make3C(const Real_IBS& _c, Kernel kernel) const   // _c is either fit face (shared Real_IBS base)
    {
        const E& aE=dynamic_cast<const E&>(*this);
        const E& cE=dynamic_cast<const E&>(_c);
        size_t Na=aE.size(), Nc=cE.size();
        PairProjection<double> s3;
        for (size_t ic=0; ic<Nc; ic++)
        {
            rsmat_t s(Na);
            for (size_t ia=0; ia<Na; ia++)
                for (size_t ib=ia; ib<Na; ib++)
                    s(ia,ib)=kernel(aE, ia, ib, cE, ic);
            s3.dense.push_back(s);
        }
        return s3;
    }
};

// --- 4-centre (HF): the ERI4 SUBSTRATE -- Direct (ab|cd) and Exchange -----------------------------
// The intricate ERI loop + symmetry packing is unchanged (hoisting it further is plan Goal D); the
// per-element integral goes through the evaluator's FourC kernel (which folds in all four norms).
// The HF CONTRACTION face (Accumulate*) is inherited from BasisSet::Orbital_ERI4_IBS -- see R1.7.
//! \brief THE UNIQUE SHELL QUARTETS of one block against itself (lane H 2026-11-11, OpenWork §4): every (ab|cd) with the 8-fold
//! permutational symmetry used ONCE, as whole shell-quartet blocks (\c FourCBlock), so the Direct (aa|cc) and the Exchange (ab|ab)
//! tensors read one set instead of evaluating 0.75 N⁴ quartets between them against N⁴/8 unique.  A value is FourC at the CANONICAL
//! permutation, so it agrees with another permutation's to roundoff (measured ≤ 1e-14 relative; tests/M_ERI4Blocks.C).  Geometry-fixed,
//! so it may be built lazily and kept (CLAUDE.md: lazy caching only for what never changes).
struct UniqueQuartets
{
    std::vector<std::pair<size_t,size_t>> S;          //!< the shells, [begin, end)
    std::vector<size_t> shellOf, offIn;               //!< per function: its shell and its place in it
    std::vector<size_t> blockOff;                     //!< per canonical quartet: its block's start in v
    std::vector<double> v;
    static size_t Pair(size_t s, size_t t) {return t*(t+1)/2+s;}   //!< s <= t
    template <class E> void Build(const E& e)
    {
        S=e.ShellRanges();
        shellOf.assign(e.size(),0); offIn.assign(e.size(),0);
        for (size_t k=0;k<S.size();k++) for (size_t f=S[k].first;f<S[k].second;f++) { shellOf[f]=k; offIn[f]=f-S[k].first; }
        const size_t nP=S.size()*(S.size()+1)/2;
        blockOff.assign(nP*(nP+1)/2, 0);
        for (size_t sb=0;sb<S.size();sb++) for (size_t sa=0;sa<=sb;sa++)
            for (size_t sd=0;sd<S.size();sd++) for (size_t sc=0;sc<=sd;sc++)
            {
                const size_t P1=Pair(sa,sb), P2=Pair(sc,sd);
                if (P2<P1) continue;
                const size_t na=S[sa].second-S[sa].first, nb=S[sb].second-S[sb].first, nc=S[sc].second-S[sc].first, nd=S[sd].second-S[sd].first;
                blockOff[P2*(P2+1)/2+P1]=v.size();
                v.resize(v.size()+na*nb*nc*nd);
                e.FourCBlock(S[sa].first,na, e,S[sb].first,nb, e,S[sc].first,nc, e,S[sd].first,nd, &v[blockOff[P2*(P2+1)/2+P1]]);
            }
    }
    double At(size_t ia, size_t ib, size_t ic, size_t id) const
    {
        size_t sa=shellOf[ia], sb=shellOf[ib], sc=shellOf[ic], sd=shellOf[id];
        size_t a=offIn[ia], b=offIn[ib], c=offIn[ic], d=offIn[id];
        if (sa>sb) { std::swap(sa,sb); std::swap(a,b); }
        if (sc>sd) { std::swap(sc,sd); std::swap(c,d); }
        size_t P1=Pair(sa,sb), P2=Pair(sc,sd);
        if (P1>P2) { std::swap(P1,P2); std::swap(sa,sc); std::swap(sb,sd); std::swap(a,c); std::swap(b,d); }
        const size_t nb=S[sb].second-S[sb].first, nc=S[sc].second-S[sc].first, nd=S[sd].second-S[sd].first;
        return v[blockOff[P2*(P2+1)/2+P1]+((a*nb+b)*nc+c)*nd+d];
    }
};

template <class E> requires (Evaluators::isHF_Evaluator<E> || Evaluators::isM_HF_Evaluator<E>)
class Orbital_ERI4_IBS
    : public virtual ::qchem::BasisSet::Orbital_ERI4_IBS<double>
    , public virtual ::qchem::BasisSet::BareCoulombSource   // the raw (m1m2|m3m4) over a subset (ACBN0, 2026-09-21)
{
public:
    // d/dR (GeometryOptimizationPlan §13.2): libcint's int2e_ip1 on the two blocks' own PGData, for the Cartesian PG_Cart evaluator when
    // the partner block is one too.  Every other pairing refuses until it has them.
    virtual bool CanPullbackERI4(const ::qchem::BasisSet::Orbital_HF_IBS<double>* other) const override
    {
        if constexpr (kLibcintERI) return dynamic_cast<const E*>(other)!=nullptr;
        else return false;
    }
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeDirectPullback(const rsmat_t& Da, const rsmat_t& Db, const ::qchem::BasisSet::Orbital_HF_IBS<double>* other,
                                                                  const Structure* st) const override {return ERI4Pullback(Da, Db, other, st, false);}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeExchangePullback(const rsmat_t& Da, const rsmat_t& Db, const ::qchem::BasisSet::Orbital_HF_IBS<double>* other,
                                                                    const Structure* st) const override {return ERI4Pullback(Da, Db, other, st, true);}
private:
    static constexpr bool kLibcintERI=std::same_as<E, Evaluators::PG_Cart_MnD::NR_Evaluator>;
    Outcome<rvec3vec_t,AD::NoGradient> ERI4Pullback(const rsmat_t& Da, const rsmat_t& Db, const ::qchem::BasisSet::Orbital_HF_IBS<double>* other,
                                                    const Structure* st, bool exchange) const
    {
        const char* what = exchange ? "exchange pullback" : "Coulomb pullback";
        if constexpr (kLibcintERI)
        {
            if (!CanPullbackERI4(other)) return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block against a partner of another family", what);
            if (!st) throw std::invalid_argument("Gaussian Orbital_ERI4_IBS: a geometry pullback needs the structure");
            using LC=Evaluators::PG_LibCint::NR_Evaluator;
            using PGD=Evaluators::PG_Cart_MnD::PGData;
            const LC a(static_cast<const PGD&>(dynamic_cast<const E&>(*this)), st), b(static_cast<const PGD&>(dynamic_cast<const E&>(*other)), st);
            auto dense=[](const rsmat_t& D){ rmat_t M(D.rows(), D.columns()); for (size_t i=0;i<D.rows();i++) for (size_t j=0;j<D.columns();j++) M(i,j)=D(i,j); return M; };
            const std::vector<rvec3_t> g=a.ERI4Gradient(b, dense(Da), dense(Db), exchange);
            rvec3vec_t F(g.size()); for (size_t A=0;A<g.size();A++) F[A]=g[A];
            return Outcome<rvec3vec_t,AD::NoGradient>::Ok(F);
        }
        else return NoGeometryDerivativeYet<rvec3vec_t>("this Gaussian block", what);
    }
public:
    //! \copydoc BasisSet::BareCoulombSource::BareCoulomb
    //! The per-element FourC kernel (all four norms folded in), the 8-fold symmetry used once per orbit.  A
    //! matrix-delivery engine (libcint) has no per-element kernel: not this increment, it throws.
    virtual ERI4Block BareCoulomb(const std::vector<size_t>& cols) const override
    {
        if constexpr (Evaluators::isM_HF_Evaluator<E>)
            throw std::runtime_error("BareCoulomb: a matrix-delivery (libcint) evaluator has no per-element (ab|cd) kernel -- not this increment");
        else
        {
            const E& e=dynamic_cast<const E&>(*this);
            const size_t m=cols.size();
            ERI4Block eri(m);
            for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++)
                for (size_t c=0;c<m;c++) for (size_t d=c;d<m;d++)
                {
                    if (a*m+b > c*m+d) continue;                     // (ab|cd) = (cd|ab)
                    const double v=e.FourC(cols[a], e, cols[b], e, cols[c], e, cols[d]);
                    eri(a,b,c,d)=eri(b,a,c,d)=eri(a,b,d,c)=eri(b,a,d,c)=v;
                    eri(c,d,a,b)=eri(d,c,a,b)=eri(c,d,b,a)=eri(d,c,b,a)=v;
                }
            return eri;
        }
    }
protected:
    // 4-centre HF Coulomb (ab|cd): a,b on this orbital basis, c,d on the partner.
    virtual ERI4 MakeDirect(const ::qchem::BasisSet::Orbital_ERI4_IBS<double>& _c) const
    {
        if constexpr (Evaluators::isM_HF_Evaluator<E>)
            return dynamic_cast<const E&>(*this).DirectMatrix(dynamic_cast<const E&>(_c));
        else if constexpr (Evaluators::isShellBlock_HF_Evaluator<E>)
        {
            const E& aE=dynamic_cast<const E&>(*this);
            const E& cE=dynamic_cast<const E&>(_c);
            if (&aE==&cE)   // THE BLOCK AGAINST ITSELF: read the unique quartets (each computed once for J and K together)
            {
                const auto held=Unique(aE); const UniqueQuartets& u=*held;
                const size_t N=aE.size();
                ERI4 J(N,N);
                for (size_t ia=0;ia<N;ia++) for (size_t ib=ia;ib<N;ib++)
                {
                    rsmat_t& Jab=J(ia,ib);
                    for (size_t ic=0;ic<N;ic++) for (size_t id=ic;id<N;id++) Jab(ic,id)=u.At(ia,ib,ic,id);
                }
                Release(/*direct*/true);
                return J;
            }
            // SHELL-BLOCKED against another block: the same entries, the same values (FourCBlock = FourC element for element).
            ERI4 J(aE.size(),cE.size());
            const auto SA=aE.ShellRanges(), SC=cE.ShellRanges();
            std::vector<double> blk;
            for (size_t sa=0;sa<SA.size();sa++) for (size_t sb=sa;sb<SA.size();sb++)
                for (size_t sc=0;sc<SC.size();sc++) for (size_t sd=sc;sd<SC.size();sd++)
                {
                    const size_t a0=SA[sa].first, na=SA[sa].second-a0, b0=SA[sb].first, nb=SA[sb].second-b0;
                    const size_t c0=SC[sc].first, nc=SC[sc].second-c0, d0=SC[sd].first, nd=SC[sd].second-d0;
                    blk.resize(na*nb*nc*nd);
                    aE.FourCBlock(a0,na, aE,b0,nb, cE,c0,nc, cE,d0,nd, blk.data());   // (a a | c c) slots
                    for (size_t a=0;a<na;a++) for (size_t b=0;b<nb;b++)
                    {
                        if (b0+b<a0+a) continue;
                        rsmat_t& Jab=J(a0+a,b0+b);
                        for (size_t c=0;c<nc;c++) for (size_t d=0;d<nd;d++)
                            if (d0+d>=c0+c) Jab(c0+c,d0+d)=blk[((a*nb+b)*nc+c)*nd+d];
                    }
                }
            return J;
        }
        else {
        const E& aE=dynamic_cast<const E&>(*this);
        const E& cE=dynamic_cast<const E&>(_c);
        size_t Na=aE.size(), Nc=cE.size();
        ERI4 J(Na,Nc);
        for (size_t ia:iv_t(0,Na))
            for (size_t ib:iv_t(ia,Na))
            {
                rsmat_t& Jab=J(ia,ib);
                for (size_t ic:iv_t(0,Nc))
                    for (size_t id:iv_t(ic,Nc))
                        Jab(ic,id)=aE.FourC(ia, aE, ib, cE, ic, cE, id);   // (a a | c c) slots
            }
        return J;
        }
    }
    // 4-centre HF Exchange: slots (a b | a b).  Symmetry packing preserved exactly.
    virtual ERI4 MakeExchange(const ::qchem::BasisSet::Orbital_ERI4_IBS<double>& _b) const
    {
        if constexpr (Evaluators::isM_HF_Evaluator<E>)
            return dynamic_cast<const E&>(*this).ExchangeMatrix(dynamic_cast<const E&>(_b));
        else if constexpr (Evaluators::isShellBlock_HF_Evaluator<E>)
        {
            const E& aE=dynamic_cast<const E&>(*this);
            const E& bE=dynamic_cast<const E&>(_b);
            if (&aE==&bE)   // THE BLOCK AGAINST ITSELF: the unique quartets, the packing below exactly
            {
                const auto held=Unique(aE); const UniqueQuartets& u=*held;
                const size_t N=aE.size();
                ERI4 K(N,N);
                for (size_t ia=0;ia<N;ia++) for (size_t ib=0;ib<N;ib++) for (size_t ic=ia;ic<N;ic++)
                {
                    rsmat_t& Kac=K(ia,ic);
                    for (size_t id=0;id<N;id++)
                    {
                        const double v=u.At(ia,ib,ic,id);
                        if (ib==id)      Kac(ib,id) = v;
                        else if (ib<id)  Kac(ib,id)+= 0.5*v;
                        else             Kac(id,ib)+= 0.5*v;
                    }
                }
                Release(/*direct*/false);
                return K;
            }
            // SHELL-BLOCKED against another block, the packing below exactly (each K entry receives at most two halves).
            ERI4 K(aE.size(),bE.size());
            const auto SA=aE.ShellRanges(), SB=bE.ShellRanges();
            std::vector<double> blk;
            for (size_t sa=0;sa<SA.size();sa++) for (size_t sc=sa;sc<SA.size();sc++)
                for (size_t sb=0;sb<SB.size();sb++) for (size_t sd=0;sd<SB.size();sd++)
                {
                    const size_t a0=SA[sa].first, na=SA[sa].second-a0, c0=SA[sc].first, nc=SA[sc].second-c0;
                    const size_t b0=SB[sb].first, nb=SB[sb].second-b0, d0=SB[sd].first, nd=SB[sd].second-d0;
                    blk.resize(na*nb*nc*nd);
                    aE.FourCBlock(a0,na, bE,b0,nb, aE,c0,nc, bE,d0,nd, blk.data());   // (a b | a b) slots
                    for (size_t a=0;a<na;a++) for (size_t c=0;c<nc;c++)
                    {
                        const size_t ia=a0+a, ic=c0+c;
                        if (ic<ia) continue;
                        rsmat_t& Kac=K(ia,ic);
                        for (size_t b=0;b<nb;b++) for (size_t d=0;d<nd;d++)
                        {
                            const size_t ib=b0+b, id=d0+d;
                            const double v=blk[((a*nb+b)*nc+c)*nd+d];
                            if (ib==id)      Kac(ib,id) = v;
                            else if (ib<id)  Kac(ib,id)+= 0.5*v;
                            else             Kac(id,ib)+= 0.5*v;
                        }
                    }
                }
            return K;
        }
        else {
        const E& aE=dynamic_cast<const E&>(*this);
        const E& bE=dynamic_cast<const E&>(_b);
        size_t Na=aE.size(), Nb=bE.size();
        ERI4 K(Na,Nb);
        for (size_t ia:iv_t(0,Na))
            for (size_t ib:iv_t(0,Nb))
                for (size_t ic:iv_t(ia,Na))
                {
                    rsmat_t& Kac=K(ia,ic);
                    for (size_t id:iv_t(0,Nb))
                    {
                        double v=aE.FourC(ia, bE, ib, aE, ic, bE, id);   // (a b | a b) slots
                        if (ib==id)      Kac(ib,id) = v;
                        else if (ib<id)  Kac(ib,id)+= 0.5*v;
                        else             Kac(id,ib)+= 0.5*v;
                    }
                }
        return K;
        }
    }
private:
    //! The unique quartets of this block against itself, built on first use under a lock (Direct and Exchange may arrive from different
    //! workers) and kept: they are a pure function of the fixed basis.
    //! A SHARED handle: the caller holds it while it reads, so \c Release (another worker's) drops only the memo's own share.
    template <class EE> std::shared_ptr<const UniqueQuartets> Unique(const EE& e) const
    {
        std::lock_guard<std::mutex> lock(itsUniqueLock->m);
        if (!itsUnique) { auto u=std::make_shared<UniqueQuartets>(); u->Build(e); itsUnique=u; itsReadJ=itsReadK=false; }
        return itsUnique;
    }
    //! Once BOTH J and K have read the quartets nothing else will (the integrals cache keeps J and K), so the memo -- about N⁴/8 doubles
    //! beside them -- is dropped; a later ask (a cache eviction) rebuilds it.
    void Release(bool direct) const
    {
        std::lock_guard<std::mutex> lock(itsUniqueLock->m);
        (direct ? itsReadJ : itsReadK)=true;
        if (itsReadJ && itsReadK) itsUnique.reset();
    }
    struct Lock { std::mutex m; };
    mutable std::shared_ptr<const UniqueQuartets> itsUnique;
    mutable bool itsReadJ=false, itsReadK=false;
    std::shared_ptr<Lock> itsUniqueLock=std::make_shared<Lock>();
};

} //namespace

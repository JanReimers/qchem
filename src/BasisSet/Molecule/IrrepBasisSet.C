// File: BasisSet/Molecule/IrrepBasisSet.C  Molecule-generic, evaluator-templated orbital IBS mixins.
//
// The molecular analog of BasisSet/Atom/IrrepBasisSet.C (module qchem.BasisSet.Atom.IBS): the i,j /
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
#include <cassert>
#include <vector>
export module qchem.BasisSet.Molecule.IBS;
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Fit_IBS;
import qchem.BasisSet.Mesh_Integrated_IBS;   // MeshIntegratorSource + MakeBeckeMeshIntegrator (field-operator)
import qchem.Mesh;                           // qcMesh::MeshParams
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Molecule.Evaluators;      // concepts + generic 1E matrix builders
export import qchem.Symmetry.Molecule.OperationRep;      // Symmetry::Molecule::AoShell (the molecule-specific 1E addition)
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

export namespace qchem::BasisSet::Molecule
{

// --- The molecular orbital 1E IBS interface -------------------------------------------------------
// The non-templated interface that ADDS molecule-specific operations to the generic Real_OIBS
// (::qchem::BasisSet::Orbital_1E_IBS<double>) -- currently just GetAoShells (the point-group SALC seam).
// A client (PG::SymmetryAdapt) can Iterate<Molecule::Orbital_1E_IBS> to reach every molecular orbital basis
// polymorphically, WITHOUT knowing its evaluator; the evaluator-templated EOrbital_1E_IBS<E> below supplies
// the evaluator-driven integral implementations.  Real-space by nature, so it lives on the molecule side and
// is absent from the dcmplx plane-wave path.
class Orbital_1E_IBS
    : public virtual ::qchem::BasisSet::Orbital_1E_IBS<double>
{
public:
    //! \brief This basis's AO shells for point-group SALC adaptation (Cartesian monomials or real solid
    //! harmonics, in the basis's own convention).  Deliveries that cannot honour a correct layout THROW
    //! (e.g. libcint-spherical, whose convention is unmatched -- S3b).
    virtual std::vector<Symmetry::Molecule::AoShell> GetAoShells() const = 0;
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
};

// --- 3-centre (DFT fit): <ab|c> for each fit component c ------------------------------------------
// Forward Cast().*ThreeC_Matrix(fit) when E delivers ERI3 (isM_DFT_Evaluator); else run the Make3C loop.
template <class E> requires (Evaluators::isDFT_Evaluator<E> || Evaluators::isM_DFT_Evaluator<E>)
class Orbital_DFT_IBS
    : public virtual ::qchem::BasisSet::Orbital_DFT_IBS<double>
    , public virtual ::qchem::BasisSet::MeshIntegratorSource<double>   // field-operator factory (CreateMeshIntegrator)
{
public:
    // Field-operator factory (the analog of CreateVxcFitBasisSet): a mesh-free Gaussian orbital basis builds
    // a fresh Becke integrator over ITSELF (the IBS IS-A VectorFunction<double>).  The static, smooth local
    // pseudopotential rides this to get <i|V_loc|j> as a raw quadrature -- no auxiliary-basis fit needed.
    virtual ::qchem::BasisSet::Mesh_Integrated_IBS<double>*
    CreateMeshIntegrator(const Structure* cl, const qcMesh::MeshParams& mp) const override
    {
        return ::qchem::BasisSet::MakeBeckeMeshIntegrator(*this, cl, mp);
    }
protected:
    virtual ERI3<double> MakeOverlap3C  (const FIT_SF_ABS& c) const
    {
        if constexpr (Evaluators::isM_DFT_Evaluator<E>)
            return dynamic_cast<const E&>(*this).OverlapThreeC_Matrix(dynamic_cast<const E&>(c));
        else return Make3C(c, [](const E& aE, size_t ia, size_t ib, const E& cE, size_t ic)
                                   {return aE.OverlapThreeC(ia, aE, ib, cE, ic);});
    }
    virtual ERI3<double> MakeRepulsion3C(const FIT_CD_ABS& c) const
    {
        if constexpr (Evaluators::isM_DFT_Evaluator<E>)
            return dynamic_cast<const E&>(*this).RepulsionThreeC_Matrix(dynamic_cast<const E&>(c));
        else return Make3C(c, [](const E& aE, size_t ia, size_t ib, const E& cE, size_t ic)
                                   {return aE.RepulsionThreeC(ia, aE, ib, cE, ic);});
    }
private:
    // For each fit component ic, build the symmetric (ia,ib) block via the supplied named 3-centre kernel
    // (which folds in all three normalizations).
    template <class Kernel>
    ERI3<double> Make3C(const Real_IBS& _c, Kernel kernel) const   // _c is either fit face (shared Real_IBS base)
    {
        const E& aE=dynamic_cast<const E&>(*this);
        const E& cE=dynamic_cast<const E&>(_c);
        size_t Na=aE.size(), Nc=cE.size();
        ERI3<double> s3;
        for (size_t ic=0; ic<Nc; ic++)
        {
            rsmat_t s(Na);
            for (size_t ia=0; ia<Na; ia++)
                for (size_t ib=ia; ib<Na; ib++)
                    s(ia,ib)=kernel(aE, ia, ib, cE, ic);
            s3.push_back(s);
        }
        return s3;
    }
};

// --- 4-centre (HF): Direct (ab|cd) and Exchange ---------------------------------------------------
// The intricate ERI loop + symmetry packing is unchanged (hoisting it further is plan Goal D); the
// per-element integral goes through the evaluator's FourC kernel (which folds in all four norms).
template <class E> requires (Evaluators::isHF_Evaluator<E> || Evaluators::isM_HF_Evaluator<E>)
class Orbital_HF_IBS
    : public virtual ::qchem::BasisSet::Orbital_HF_IBS<double>
{
protected:
    // 4-centre HF Coulomb (ab|cd): a,b on this orbital basis, c,d on the partner.
    virtual ERI4 MakeDirect(const ::qchem::BasisSet::Orbital_HF_IBS<double>& _c) const
    {
        if constexpr (Evaluators::isM_HF_Evaluator<E>)
            return dynamic_cast<const E&>(*this).DirectMatrix(dynamic_cast<const E&>(_c));
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
    virtual ERI4 MakeExchange(const ::qchem::BasisSet::Orbital_HF_IBS<double>& _b) const
    {
        if constexpr (Evaluators::isM_HF_Evaluator<E>)
            return dynamic_cast<const E&>(*this).ExchangeMatrix(dynamic_cast<const E&>(_b));
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
};

} //namespace

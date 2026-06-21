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
// bodies are exactly the ones that used to be hand-written in PolarizedGaussian/Imp/Orbital_IBS.C.
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
import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Molecule.Evaluators;      // concepts + generic 1E matrix builders
import qchem.Cluster;
import qchem.Types;
import qchem.Blaze;

export namespace BasisSet::Molecule
{

// --- 1E: Overlap / Kinetic(<p^2>) / Nuclear -------------------------------------------------------
// The element kernels and the i,j loop both live on the evaluator side (Evaluators::*Matrix); this mixin
// just forwards the IBS's MakeXxx() virtuals onto them.
template <Evaluators::is1E_Evaluator E> class Orbital_1E_IBS
    : public virtual ::BasisSet::Orbital_1E_IBS<double>
{
protected:
    virtual rsmat_t MakeOverlap(                 ) const {return Evaluators::OverlapMatrix(Cast()    );}
    virtual rsmat_t MakeKinetic(                 ) const {return Evaluators::KineticMatrix(Cast()    );}
    virtual rsmat_t MakeNuclear(const Cluster* cl) const {return Evaluators::NuclearMatrix(Cast(), cl);}
    const E& Cast() const {return dynamic_cast<const E&>(*this);}
};

// --- 3-centre (DFT fit): <ab|c> for each fit component c ------------------------------------------
template <Evaluators::isDFT_Evaluator E> class Orbital_DFT_IBS
    : public virtual ::BasisSet::Orbital_DFT_IBS<double>
{
protected:
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const
    {
        return Make3C(c, [](const E& aE, size_t ia, size_t ib, const E& cE, size_t ic)
                            {return aE.OverlapThreeC(ia, aE, ib, cE, ic);});
    }
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const
    {
        return Make3C(c, [](const E& aE, size_t ia, size_t ib, const E& cE, size_t ic)
                            {return aE.RepulsionThreeC(ia, aE, ib, cE, ic);});
    }
private:
    // For each fit component ic, build the symmetric (ia,ib) block via the supplied named 3-centre kernel
    // (which folds in all three normalizations).
    template <class Kernel>
    ERI3<double> Make3C(const Fit_IBS& _c, Kernel kernel) const
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
template <Evaluators::isHF_Evaluator E> class Orbital_HF_IBS
    : public virtual ::BasisSet::Orbital_HF_IBS<double>
{
protected:
    // 4-centre HF Coulomb (ab|cd): a,b on this orbital basis, c,d on the partner.
    virtual ERI4 MakeDirect(const ::BasisSet::Orbital_HF_IBS<double>& _c) const
    {
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
    // 4-centre HF Exchange: slots (a b | a b).  Symmetry packing preserved exactly.
    virtual ERI4 MakeExchange(const ::BasisSet::Orbital_HF_IBS<double>& _b) const
    {
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
};

} //namespace

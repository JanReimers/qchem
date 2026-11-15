// File: BasisSet/Radial/IBS.C Atom specific irrep basis sets.
module;
#include <iosfwd>
#include <memory>
#include <cassert>
#include <algorithm>
#include <iostream>
#include "forward.H"

export module qchem.BasisSet.Radial.IBS;
export import qchem.BasisSet.BasisCache;   // a block's constructor takes its parent's cache (IntegralCacheOwnershipPlan S3a)
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.Internal.Orbital_DHF_IBS;
import qchem.BasisSet.IrrepBasisSet;
export import qchem.BasisSet.ImplicitAngular_IBS;   // the radial/implicit-Y_lm capability atoms realize
import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.Fit_IBS;   // the shared Gaussian-family implementation base we DERIVE
import qchem.BasisSet.Internal.Orbital_ERI4_IBS;
import qchem.BasisSet.Internal.Cache4;   // Cache4 / Cacheable4: the parent's Rk table (S5a)
import qchem.BasisSet.Radial.Evaluators;
import qchem.BasisSet.Internal.Cache4;
import qchem.Symmetry.Factory;
import qchem.Constants;
import qchem.Blaze;
import qchem.Reporting;   // a truncated Coulomb metric reports itself

export namespace qchem::BasisSet
{
namespace Radial
{

    using namespace Evaluators;
    using EvaluatorsBase = Evaluators::Evaluator;

//  Atom-only identity face.  An atom's cache identity splits into a RADIAL part (exponents / contraction /
//  grid) and an ANGULAR part (the l,m channel); the radial part alone also keys the Slater Rk cache
//  (Cache4), which is why atoms keep the two separable.  Centre is pinned at the nucleus, so radial|angular
//  is a complete identity.  This lives in BasisSet/Atom because it is atom-specific: the structure-neutral
//  BasisSet layer knows only the single BasisSetID (a neutral class that reaches for RadialID/AngularID has
//  a design flaw -- it should compose from BasisSetID instead; cf. the RKB Orbital_RKB_IBS_Imp).
class RadialAngularID : public virtual BasisSet::IrrepBasisSet_IDs
{
public:
    virtual std::string  RadialID() const=0;
    virtual std::string AngularID() const=0;
    virtual std::string BasisSetID() const override {return RadialID()+"|"+AngularID();}
};

//
//  Common IrrepBasisSet functionality for atom basis sets.  All the work is done by the evaluator
//
// ⚠ THE SURVIVING FAKE (doc/CleanupCandidates.md R1.0).  This block declares itself Evaluatable_IBS, and
// its op(r) returns the purely RADIAL chi_i(|r|) with the irrep's Y_lm silently omitted -- so the promise
// is kept in form and broken in substance.  It is CONTAINED, not cured: ImplicitAngular_IBS makes the
// fakeness explicit in the type system and PP_NonLocal cross-casts to it for a per-l radial assembly.
// The cure is R1.0 step (1)'s remaining conversions -- PP_Local's MatrixOverlap and PP_NonLocal's
// explicit-angular Overlap become MakeOverlap/MakeOverlap3C on the basis, and the atom density/orbital
// paths follow -- after which this derivation drops and nothing evaluates an atom block from outside.
// (Step (2), an HONEST op(r) with real Y_lm, is explicitly out of scope -- user 2026-08-22.)
template <isOpr_Evaluator E> class IrrepBasisSetImp
    : public virtual BasisSet::Evaluatable_IBS<double>
    , public virtual RadialAngularID   // atom identity: RadialID/AngularID compose BasisSetID
    , public virtual BasisSet::ImplicitAngular_IBS  // the Y_lm is carried by the IRREP, not the function
    , public BasisSet::IrrepBasisSetImp<double> //Pulls in Symmetry support
    , public Symmetry::IdentityRepresentation<double>   // an l block stores the RADIAL factor: identity under every op fixing
                                                         // the nucleus (D^l acts on the implicit m partners) -- V-BS-SYMMETRY S3
{
public:
    IrrepBasisSetImp(const std::shared_ptr<BasisSet::BasisCache>& cache, const sym_t& yl) : BasisSet::IrrepBasisSetImp<double>(yl, cache)
        , Symmetry::IdentityRepresentation<double>(Symmetry::Stabiliser::FixingPoint(rvec3_t(0,0,0))) {}

    virtual size_t GetNumFunctions() const override {return Cast().size();}
    // using statements in the final class don't seem to work, so we need to function forward.
    virtual rvec_t     operator() (const rvec3_t& r) const override {return Cast().operator()(r);}
    // d/dR, forward (GeometryOptimizationPlan §13.2): the gradient source awaits the user's ruling, so every family refuses for now.
    virtual bool CanGiveJets() const override {return false;}
    virtual Outcome<size_t,AD::NoGradient> JetsOn(const rvec3vec_t&, size_t, const Structure*, BasisSet::PointJets<double>&) const override
        {return BasisSet::NoGeometryDerivativeYet<size_t>("a radial (atom) block", "function jets");}
    virtual rvec3vec_t Gradient   (const rvec3_t& r) const override {return Cast().Gradient  (r);}

    // ImplicitAngular_IBS -- an atom block IS radial: its functions omit the irrep's Y_lm (see the face).
    virtual int    ImplicitL   ()         const override {return Cast().Getl();}
    virtual rvec_t RadialValues(double r) const override {return Cast().operator()(rvec3_t(r,0,0));}
   
    virtual std::string RadialID () const override {return Cast().RadialID();}
    virtual std::string AngularID() const override {return Cast().AngularID();}
    virtual std::string Name     () const override {return Cast().Name();}
protected:
    auto& Cast() const {return dynamic_cast<const E&>(*this);}
};

//
//  Implement these integral engines separately as they shared between NR and RKB 1E orbital IBS implementations.
//  Fit_IBS also uses Overlap.
//
template <is1E_Evaluator E> class Integrals_Overlap
: public virtual BasisSet::Integrals_Overlap<double>
{
public:
    // d/dR (GeometryOptimizationPlan §13.2): no geometry derivative in this family yet.
    virtual bool CanPullbackOverlap() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<double>&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a radial (atom) block", "overlap pullback");}
protected:
    virtual smat_t<double> MakeOverlap() const 
    {
        auto& e=dynamic_cast<const E&>(*this);
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Overlap(i,j);

        return S;
    }
};
template <is1E_Evaluator E> class Integrals_Kinetic
: public virtual BasisSet::Integrals_Kinetic<double>
{
public:
    // d/dR (GeometryOptimizationPlan §13.2): no geometry derivative in this family yet.
    virtual bool CanPullbackKinetic() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeKineticPullback(const hmat_t<double>&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a radial (atom) block", "kinetic pullback");}
protected:
    // Builds the kinetic BUILDING BLOCK \f$\langle p^2\rangle=\langle -\nabla^2\rangle\f$ (NOT energy,
    // no 1/2 -- see BasisSet/Orbital_1E_IBS.C).  For atoms \f$-\nabla^2\f$ splits into the radial
    // gradient piece e.Grad2(i,j) plus the centrifugal (angular) term \f$l(l+1)\langle r^{-2}\rangle\f$.
    virtual smat_t<double> MakeKinetic() const
    {
        auto& e=dynamic_cast<const E&>(*this);
        int l=e.Getl();
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Grad2(i,j) + l*(l+1)*e.Inv_r2(i,j);

        return S;
    }
};
template <is1E_Evaluator E> class Integrals_Nuclear
: public virtual BasisSet::Integrals_Nuclear<double>
{
public:
    // d/dR (GeometryOptimizationPlan §13.2): no geometry derivative in this family yet.
    virtual bool CanPullbackNuclear() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeNuclearPullback(const hmat_t<double>&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a radial (atom) block", "nuclear pullback");}
protected:
    virtual smat_t<double> MakeNuclear(const Structure* cl) const 
    {
        assert(cl);
        assert(cl->GetNumAtoms()==1); //This supposed to be an atom after all!
        int Z=-cl->GetNuclearCharge(); 
        auto& e=dynamic_cast<const E&>(*this);
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= Z*e.Inv_r1(i,j);

        return S;
    }
};

// Use E prefix to avoid name clash with the interface class Fit_IBS
template <isFit_Evaluator Evaluator> class EFit_IBS
    : public virtual Fit_IBS 
    , public Integrals_Overlap<Evaluator>
    , public IrrepBasisSetImp<Evaluator>
    , public Evaluator
{
    using IrrepBasisSetImp<Evaluator>::Cast;
public:
    //! The Coulomb fit: analytic, no mesh.
    EFit_IBS(const std::shared_ptr<BasisSet::BasisCache>& cache, const Evaluator& e)
        : Fit_IBS(), IrrepBasisSetImp<Evaluator>(cache, Symmetry::YFactory()), Evaluator(e) {};
    //! The Vxc fit: its projection mesh, built from \a st and \a mp.
    EFit_IBS(const std::shared_ptr<BasisSet::BasisCache>& cache, const Evaluator& e, const Structure& st, const qcMesh::MeshParams& mp)
        : Fit_IBS(st,mp), IrrepBasisSetImp<Evaluator>(cache, Symmetry::YFactory()), Evaluator(e) {};
    // d/dR (GeometryOptimizationPlan §13.2): a one-centre (atom) fit has no geometry to move.
    virtual bool CanPullbackRepulsion() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const rsmat_t&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("the radial fit block", "Coulomb-metric pullback");}

    //! \brief The Coulomb metric's inverse, ROBUST to a near-singular metric.  A dense Slater s set (the atom pilot's raised Coulomb-fit rungs,
    //! \c CoulombFitLadder) makes G near-singular, and a plain inverse then hands the Dunlap fit garbage (U, rung 2: |ΔE| 7e4 Ha, lane B
    //! 2026-11-14).  So G = V diag(w) V^T, and the modes with w < \c kMetricRelTol x w_max are DROPPED: the fit then lives on the span G
    //! resolves, and Dunlap's robust energy, stationary in the coefficients, stays second order in the residual on that span.  A metric
    //! with no such mode takes the plain inverse, unchanged to the bit.  MEASURED 2026-11-14 (Kr-U, High pools): it is not only the raised
    //! rungs -- the classic fit (rung 0) is near-singular too, 4-9 of 36-46 modes on Gaussian and 1 of 21-22 on Slater Au, Rn, U, which
    //! the plain inverse had been inverting silently; the dense Slater rungs drop most of theirs (U rung 4: 472 of 510).  Truncated, every
    //! Dunlap energy agrees with the radial-Poisson Hartree to a few μHa (A_DFT_heavy).  A truncation REPORTS itself (count and floor).
    static constexpr double kMetricRelTol=1e-12;
    virtual rsmat_t MakeInvRepulsion() const override
    {
        const rsmat_t G=MakeRepulsion();
        const size_t n=G.rows();
        rvec_t w; rmat_t V;
        blazem::eigen(G, w, V);
        double wmax=0.0;
        for (size_t k=0;k<n;k++) wmax=std::max(wmax, w[k]);
        bool drop=false;
        for (size_t k=0;k<n;k++) drop = drop || w[k]<kMetricRelTol*wmax;
        if (!drop) return blazem::inv(G);
        size_t dropped=0;
        for (size_t k=0;k<n;k++) if (w[k]<kMetricRelTol*wmax) dropped++;
        // A truncation is a choice the run must show (a silent one is a substitution): the count and the floor, one level in, under
        // whatever scope the caller opened (the Hartree fit's), and on the console.
        qchem::report::Set("coulombMetricModesDropped", (long)dropped);
        qchem::report::Set("coulombMetricFloor", kMetricRelTol*wmax);
        std::cout<<"[coulomb fit] the metric of "<<n<<" functions is near-singular: "<<dropped<<" mode(s) below the floor "
                 <<kMetricRelTol*wmax<<" (= "<<kMetricRelTol<<" x the largest eigenvalue) are dropped"<<std::endl;
        rsmat_t Ginv(n);
        for (size_t i=0;i<n;i++)
            for (size_t j=i;j<n;j++)
            {
                double s=0.0;
                for (size_t k=0;k<n;k++) if (w[k]>=kMetricRelTol*wmax) s+=V(i,k)*V(j,k)/w[k];
                Ginv(i,j)=s;
            }
        return Ginv;
    }
    virtual rsmat_t MakeRepulsion(                ) const 
    {
        auto& e=Cast();
        rsmat_t S(e.size());
        for (auto i:e.indices())
            for (auto j:e.indices(i))
                S(i,j)= e.Repulsion(i,j);

        return S;
    }
    virtual  rmat_t MakeRepulsion(const rFIT_CD_ABS& f) const
    {
        auto& ea=Cast();
        auto& eb=dynamic_cast<const Evaluator&>(f);
        rmat_t S(ea.size(),eb.size());
        for (auto i:ea.indices())
            for (auto j:eb.indices())
                S(i,j)= ea.Repulsion(i,j,eb);

        return S;
    }
    virtual  rvec_t MakeCharge   (                ) const 
    {
        auto& e=Cast();
        rvec_t c(e.size());
        for (auto i:e.indices())
            c[i]=e.Charge(i);
        return c;
    }
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << "Atom fit IBS " << Name() << " " << GetSymmetry();
    }

};

//
//  1E NR orbital for atoms.  Use mixins to get the integral evaluations.
//
template <is1E_Evaluator E> class Orbital_1E_IBS
    : public virtual BasisSet::Orbital_1E_IBS<double> //This part has the symmetry.
    , public Integrals_Overlap<E>
    , public Integrals_Kinetic<E>
    , public Integrals_Nuclear<E>
{
public:
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << Name() << " " << GetSymmetry();
    }
};


template <isDFT_Evaluator E> class Orbital_DFT_IBS
    : public virtual BasisSet::Orbital_DFT_IBS<double>
{
protected:
    virtual PairProjection<double> MakeOverlap3C  (const rFIT_SF_ABS& _c) const
    {
        auto& ab=dynamic_cast<const E&>(*this);
        auto& c =dynamic_cast<const E&>(_c);
        PairProjection<double> S3;
        for (auto ic:c.indices())
        {
            rsmat_t S(ab.size());
            for (auto i:ab.indices())
                for (auto j:ab.indices(i))
                    S(i,j)=ab.Overlap(i,j,c,ic);
            
            S3.dense.push_back(S);
        }
        return S3;

    }
public:
    // d/dR (GeometryOptimizationPlan §13.2): a one-centre (atom) block has no geometry to move.
    virtual bool CanPullbackRepulsion3C(const rFIT_CD_ABS&) const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsion3CPullback(const rFIT_CD_ABS&, const hmat_t<double>&, const rvec_t&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a radial (atom) block", "3-centre repulsion pullback");}
protected:
    virtual PairProjection<double> MakeRepulsion3C(const rFIT_CD_ABS& _c) const
    {
        auto& ab=dynamic_cast<const E&>(*this);
        auto& c =dynamic_cast<const E&>(_c);
        PairProjection<double> S3;
        for (auto ic:c.indices()) 
        {
            rsmat_t S(ab.size());
            for (auto i:ab.indices())
                for (auto j:ab.indices(i))
                    S(i,j)=ab.Repulsion(i,j,c,ic);  
            
            S3.dense.push_back(S);
        }
        return S3;

    }
};




// The ERI4 SUBSTRATE for an atomic irrep block: the Cache4/Grouper/Ak radial-Slater design.  The HF
// CONTRACTION face (Accumulate*) is inherited from BasisSet::Orbital_ERI4_IBS -- see R1.7.
template <isHF_Evaluator E> class Orbital_ERI4_IBS
    : public virtual BasisSet::Orbital_ERI4_IBS<double> 

{
public:
    // d/dR (GeometryOptimizationPlan §13.2): a one-centre (atom) block has no geometry to move.
    virtual bool CanPullbackERI4(const BasisSet::Orbital_HF_IBS<double>*) const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeDirectPullback(const rsmat_t&, const rsmat_t&, const BasisSet::Orbital_HF_IBS<double>*, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a radial (atom) block", "Coulomb pullback");}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeExchangePullback(const rsmat_t&, const rsmat_t&, const BasisSet::Orbital_HF_IBS<double>*, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("a radial (atom) block", "exchange pullback");}
protected:
    virtual ERI4 MakeDirect  (const BasisSet::Orbital_ERI4_IBS<double>& _c) const;
    virtual ERI4 MakeExchange(const BasisSet::Orbital_ERI4_IBS<double>& _c) const;

};

template <is1E_Evaluator E> class Orbital_RKBL_IBS
    : public virtual BasisSet::Orbital_RKBL_IBS<double> 
    , public Integrals_Overlap<E>
    , public Integrals_Kinetic<E> //Used internally only.
    , public Integrals_Nuclear<E>
{
public:
    // RELATIVISTIC kinetic (RKB L/S cross term).  Integrals_Kinetic<E>::MakeKinetic() is the
    // \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block (no 1/2); here it is scaled by 1/(2 c_light).
    // WARNING: this 1/(2 c_light) and the c_light in Imp/Orbital_DHF_IBS.C::MakeKinetic() multiply out
    // to ~1/2 -- the c's cancel.  Whether each factor sits in the *right* place is UNVERIFIED (see the
    // \warning in BasisSet/Internal/Orbital_DHF_IBS.C); kept as-is, guarded by A_DHF.
    virtual rmat_t  MakeKinetic(const Orbital_RKBS_IBS<double>& rkbs) const
    {
        auto& ea=dynamic_cast<const E&>(*this);
        auto& eb=dynamic_cast<const E&>(rkbs); //The RKBS evaluator should be derived from the RKBL evaluator.
        assert(ea.Getl()==eb.Getl());
        assert(ea.size()==eb.size());
        return Integrals_Kinetic<E>::MakeKinetic()/(2*c_light);
    }
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << Name() << " " << GetSymmetry();
    }
};

template <is1E_Evaluator E> class Orbital_RKBS_IBS
    : public virtual BasisSet::Orbital_RKBS_IBS<double> 
    , public Integrals_Kinetic<E>
    , public Integrals_Nuclear<E> //RKBS Evaluator overrides Inv_r1 definition
{
public:
    // RKB small-component OVERLAP.  In RKB the small-component metric \f$\langle\chi^S|\chi^S\rangle\f$
    // is proportional to \f$\langle p^2\rangle\f$, so this reuses the kinetic \f$\langle p^2\rangle\f$
    // block (Integrals_Kinetic::MakeKinetic, no 1/2) as the overlap.  Any \f$1/(4c^2)\f$ normally
    // expected here is absent -- presumed to cancel against RestMass/scaling elsewhere; UNVERIFIED
    // (see the \warning in BasisSet/Internal/Orbital_DHF_IBS.C).
public:
    // d/dR (GeometryOptimizationPlan §13.2): no geometry derivative in this family yet (its MakeOverlap is the small-component metric).
    virtual bool CanPullbackOverlap() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<double>&, const Structure*) const override
        {return BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("the radial RKB small-component block", "overlap pullback");}
    virtual rsmat_t MakeOverlap() const
    {
        return Integrals_Kinetic<E>::MakeKinetic();
    }
    virtual std::ostream&  Write(std::ostream& os) const
    {
        return os << Name() << " " << GetSymmetry();
    }

};



template <isHF_Evaluator E> ERI4 Orbital_ERI4_IBS<E>::MakeDirect(const BasisSet::Orbital_ERI4_IBS<double>& _c) const 
{
    auto& a=dynamic_cast<const E&>(*this);
    auto& c=dynamic_cast<const E&>(_c);
    int la=a.Getl(), lc=c.Getl();
    assert(a.RadialType()==c.RadialType());
    size_t spanab=a.maxSpan(),spancd=c.maxSpan();
    size_t Na=a.size(), Nc=c.size();
    rvec11_t Akac=a.DirectAk(c);
    const Cache4* Rk_cache=&this->Cache().GetCache4(a.RadialType());   // the parent's (S5a); c is a block of the same parent
    ERI4 J(Na,Nc);

    for (size_t ia:a.indices())
    {
        Rk_cache->loop_1(a.es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c.indices())
        {
            Rk_cache->loop_2(c.es_index(ic));
            for (size_t ib:a.indices())
            {
                if (ib<ia) continue; 
                if (ib>ia+spanab) continue;
                rsmat_t& Jab=J(ia,ib);
                Rk_cache->loop_3(a.es_index(ib));
                for (size_t id:c.indices())
                {
                    if (id<ic) continue;
                    if (id>ic+spancd) continue;
                    if (Jab(ic,id)!=0.0)
                    {
                        std::cout << "overwriting Jnew(" << ia << " " << ib << " " << ic << " " << id << ")="; 
                        std::cout << Jab(ic,id) << std::endl;    
                        assert(false);
                    }
                    double norm=a.Norm(ia)*a.Norm(ib)*c.Norm(ic)*c.Norm(id);
                    const Cacheable4* rk=Rk_cache->loop_4(c.es_index(id));
                    double AkacRkac=E::direct(rk,la,lc,Akac);  //static call.
                    Jab(ic,id)=AkacRkac*norm;
                }
            }
        }
    }
    return J;
}

template <isHF_Evaluator E> ERI4 Orbital_ERI4_IBS<E>::MakeExchange(const BasisSet::Orbital_ERI4_IBS<double>& _c) const 
{
    auto& a=dynamic_cast<const E&>(*this);
    auto& c=dynamic_cast<const E&>(_c);
    assert(a.RadialType()==c.RadialType());
    size_t spanab=a.maxSpan(),spancd=c.maxSpan();
    size_t Na=a.size(), Nc=c.size();
    int la=a.Getl(), lc=c.Getl();
    rvec11_t Akac=a.ExchangeAk(c);
    const Cache4* Rk_cache=&this->Cache().GetCache4(a.RadialType());   // the parent's (S5a); c is a block of the same parent

    ERI4 K(Na,Nc);
    for (size_t ia:a.indices())
    {
        Rk_cache->loop_1(a.es_index(ia)); //Start a cache for Gaussian::RkEngine*
        for (size_t ic:c.indices())
        {
            if (ia>ic+spanab || ic>ia+spancd) continue;
            for (size_t ib:a.indices(ia))
            {
                rsmat_t& Kab=K(ia,ib);
                Rk_cache->loop_2(a.es_index(ib));
                Rk_cache->loop_3(c.es_index(ic));
                for (size_t id:c.indices())
                {
                    if (id>ib+spancd || ib>id+spanab) continue;
                    double norm=a.Norm(ia)*a.Norm(ib)*c.Norm(ic)*c.Norm(id);
                    const Cacheable4* rk=Rk_cache->loop_4(c.es_index(id));
                    double AkacRKac=E::exchange(rk,la,lc,Akac); //static call.
                    if (ic==id)
                        Kab(ic,id)=AkacRKac*norm; 
                    else if (id<ic)
                        Kab(id,ic)+=0.5*AkacRKac*norm; 
                    else
                        Kab(ic,id)+=0.5*AkacRKac*norm; 

                }
            }
        }
    }

    return K;

}


}} //namespaces
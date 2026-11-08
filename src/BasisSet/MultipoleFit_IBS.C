// File: BasisSet/MultipoleFit_IBS.C  The multipole fit basis: radial channels per (site, l, m) on an atom-centred mesh.
//
// f(r) ~ sum_A sum_lm f^A_lm(|r-R_A|) Y_lm, with f^A_lm(r_i) = int dOmega Y_lm p_A f (Becke & Dickson 1988; DMol; FHI-aims).  The
// functions are one per (site, shell, l, m), and the expansion is a real function everywhere: a spline in r times Y_lm.
//
// A DECORATOR OF THE δ BASIS.  Both ends of every contraction are the mesh's, so this class puts a transform on each side of the δ
// basis's own forward and adjoint and adds no table of its own:
//     forward   D -> rho(r_g) -> moments / shell weight        the channels
//     adjoint   channels -> f(r_g) -> sum_g w_g f_g chi_i chi_j
// The projection needs no solve, since the partition of unity gives sum_A p_A f = f before truncation, so isOrtho() is true in the
// sense it is used: no metric solve.  It is a dual pair, not an orthogonal set.  The projection carries p_A and the expansion does
// not, and because a site's functions overlap its neighbours' the metric is not literally diagonal.  OverlapDiagonal is the shell
// weight, the divisor that turns a projection into a coefficient.
//
// IT IS A FINITE-STRUCTURE BASIS.  Hamiltonian::AdmitVxc and AdmitCoulomb refuse it on a periodic one.
module;
#include <cassert>
#include <cmath>
#include <map>
#include <memory>
#include <mutex>
#include <numbers>
#include <ostream>
#include <string>
export module qchem.BasisSet.MultipoleFit_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.Mesh.SiteMultipoles;
export import qchem.BasisSet.CoulombDeltaFit_IBS;   // PointCoulombIntegrals: this basis answers its own Coulomb integrals
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Blaze;
import qchem.Parallel;   // LazyLock -- the per-block cache's guard
import qchem.Symmetry;
import qchem.Types;

export namespace qchem::BasisSet
{

//! \brief The multipole realization of the scalar fit face, over a \f$\delta\f$ basis on the same mesh.
class MultipoleFit_IBS
    : public virtual rFIT_SF_ABS
    , public virtual PointCoulombIntegrals   //!< its own radial-Poisson solve over the sites (V-ONE-OF-N 3b-ii, ruling (f))
    , public         IrrepBasisSetImp<double>
    , public         Symmetry::IdentityRepresentation<double>
{
public:
    //! \a delta is the \f$\delta\f$ basis on the mesh \a sites was built from, and \a sym the fit lineage's irrep.
    MultipoleFit_IBS(std::shared_ptr<const rFIT_SF_ABS> delta, std::shared_ptr<const qcMesh::SiteMultipoles> sites, const sym_t& sym)
        : IrrepBasisSetImp<double>(sym), itsDelta(std::move(delta)), itsSites(std::move(sites))
        , itsW(itsDelta->Charge()), itsDiagonal(MakeDiagonal())
    {}
    // --- its Coulomb integrals: the multipole solve over its own sites, with its own metric (the shell weights)
    rvec_t Potential    (const rvec_t& q) const override                               {return Coulomb().Potential(q);}
    double SelfRepulsion(const rvec_t& q, const rvec_t& v) const override              {return Coulomb().SelfRepulsion(q,v);}
    double HartreeEnergy(double trDV, const rvec_t& q, const rvec_t& v) const override {return Coulomb().HartreeEnergy(trDV,q,v);}
    bool   UsesTrDV     () const override                                              {return Coulomb().UsesTrDV();}

    bool isOrtho() const override {return true;}   //!< no metric solve: see the file header

    //! \f$\langle f_a|f\rangle\f$: the moments of the \f$\delta\f$ projection, which is \f$w_gf(r_g)\f$.
    rvec_t Overlap(const ScalarFunction<double>& f) const override {return itsSites->Flatten(itsSites->Moments(itsDelta->Overlap(f)));}
    //! The shell weight \f$r_i^2dr_i\f$: \f$\langle f_a|f\rangle\f$ over this is \f$f^A_{lm}(r_i)\f$.
    rvec_t OverlapDiagonal() const override {return itsDiagonal;}
    //! Only \f$l=0\f$ carries charge, \f$\sqrt{4\pi}\f$ times its shell weight.  So the fit keeps the mesh's charge exactly.
    rvec_t Charge() const override
    {
        rvec_t c(itsDiagonal.size(),0.0);
        const size_t nLM=itsSites->NumLM();
        for (size_t a=0; a<c.size(); a+=nLM) c[a]=2.0*std::sqrt(std::numbers::pi)*itsDiagonal[a];
        return c;
    }
    //! \copydoc BasisSet::Integrals_Overlap3C::Overlap3C
    const PairProjection<double>& Overlap3C(const Orbital_DFT_IBS<double,double>& orb) const override
    {
        const PairProjection<double>& d=itsDelta->Overlap3C(orb);   // the δ basis's own lock, taken before mine
        std::lock_guard<std::mutex> lock(itsLock.m);
        const sym_t& id=orb.GetSymt();
        auto it=itsO3.find(id);
        if (it!=itsO3.end()) return it->second;
        PairProjection<double> g;
        g.applyRaw       =[this,&d](const hmat_t<double>& D) {return Channels(d.applyRaw(D));};
        g.applyRawAdjoint=[this,&d](const rvec_t& c) {return d.applyRawAdjoint(itsSites->Synthesize(itsSites->Split(c)));};
        if (d.applyRawFactored) g.applyRawFactored=[this,&d](const mat_t<double>& L) {return Channels(d.applyRawFactored(L));};
        return itsO3.emplace(id,std::move(g)).first->second;
    }

    size_t GetNumFunctions() const override {return itsSites->NumChannels();}
    std::string   Name      () const override {return "MultipoleFit";}
    std::string   BasisSetID() const override
        {return Name()+"|channels="+std::to_string(GetNumFunctions())+"|lMax="+std::to_string(itsSites->LMax());}
    std::ostream& Write     (std::ostream& os) const override
        {return os << Name() << " fit IBS: " << GetNumFunctions() << " channels to l=" << itsSites->LMax();}

private:
    rvec_t MakeDiagonal() const
    {
        rvec_t d(itsSites->NumChannels());
        const auto& W=itsSites->ShellWeights();
        const size_t nLM=itsSites->NumLM();
        size_t k=0;
        for (size_t a=0; a<itsSites->NumSites(); a++)
            for (size_t i=0; i<W.size(); i++)
                for (size_t c=0; c<nLM; c++) d[k++]=W[i];
        return d;
    }
    //! A density's values at the points, as its coefficients on my functions.
    //! The solve, built on first use: a pure function of the sites and the metric (fixed), and an XC-only fit basis never needs it.
    const qcMesh::CoulombOperator& Coulomb() const
    {
        std::lock_guard<std::mutex> lock(itsLock.m);
        if (!itsCoulomb) itsCoulomb=qcMesh::MakeMultipoleCoulomb(itsSites, itsDiagonal);
        return *itsCoulomb;
    }
    rvec_t Channels(const rvec_t& rho) const
    {
        rvec_t q(rho.size());
        for (size_t g=0; g<q.size(); g++) q[g]=itsW[g]*rho[g];
        rvec_t c=itsSites->Flatten(itsSites->Moments(q));
        for (size_t a=0; a<c.size(); a++) c[a]/=itsDiagonal[a];
        return c;
    }

    std::shared_ptr<const rFIT_SF_ABS>            itsDelta;
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    rvec_t itsW;          //!< the mesh's weights
    rvec_t itsDiagonal;   //!< the shell weight of each channel
    mutable std::unique_ptr<qcMesh::CoulombOperator> itsCoulomb;   //!< its Coulomb integrals, built on first use (Coulomb())
    mutable SymMap<PairProjection<double>> itsO3;   //!< the composed pairs, by SPATIAL symmetry
    qchem::LazyLock itsLock;
};

} //namespace

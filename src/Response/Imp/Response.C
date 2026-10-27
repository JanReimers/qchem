// File: Response/Imp/Response.C  Walk a converged Bloch wave function into the response vocabulary.
module;
#include <iomanip>
#include <cmath>
#include <complex>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
module qchem.Response;
import qchem.Orbitals;                   // TOrbitals / TOrbital (the per-block orbital sets)
import qchem.BasisSet.Orbital_DFT_IBS;   // the block basis the projector amplitudes are asked on
import qchem.Mesh.Quadrature;            // MatrixOverlap: <chi_a| r_i |chi_b> on a mesh (the dipole probe)
import qchem.ScalarFunction;
import qchem.Blaze;

namespace qchem::Response
{

namespace {

// Every block is visited in the wave function's own QN order, by both adapters -- so block b of the
// Reference and block b of the probe are the same (k, σ).  A block's orbitals are either scalar (a real
// TRIM block inside a complex run is TOrbitals<double>); the visitor hands each branch its native type.
template <class T, class Visit> void ForEachBlock(const WaveFunction::tWaveFunction<T>& wf, Visit&& visit)
{
    for (const Irrep& ir : wf.GetQNs())
    {
        const Orbitals::Orbitals* os=wf.GetOrbitals(ir);
        if (const auto* c=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(os)) visit(ir, *c);
        else if (const auto* r=dynamic_cast<const Orbitals::TOrbitals<double>*>(os)) visit(ir, *r);
        else throw std::logic_error("Response: an orbital block of unknown scalar");
    }
}

template <class U> mat_t<U> Coefficients(const Orbitals::TOrbitals<U>& os, size_t nbasis)
{
    const size_t norb=os.GetNumOrbitals();
    mat_t<U> C(nbasis, norb);
    size_t i=0;
    for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>())
    {
        const vec_t<U>& c=o->GetCoeff();
        for (size_t a=0;a<c.size();a++) C(a,i)=c[a];
        i++;
    }
    return C;
}

cmat_t ToComplex(const mat_t<double>& m) {cmat_t c(m.rows(),m.columns()); for (size_t i=0;i<m.rows();i++) for (size_t j=0;j<m.columns();j++) c(i,j)=m(i,j); return c;}
cmat_t ToComplex(const mat_t<dcmplx>& m) {return m;}

//! One Cartesian coordinate \f$r_i\f$ as a scalar field: the dipole operator's multiplicative kernel.
class Coordinate : public ScalarFunction<double>
{
public:
    explicit Coordinate(int i) : itsI(i) {}
    virtual double  operator()(const rvec3_t& r) const override {return itsI==0 ? r.x : itsI==1 ? r.y : r.z;}
    virtual rvec3_t Gradient  (const rvec3_t&  ) const override {return rvec3_t(itsI==0, itsI==1, itsI==2);}
private:
    int itsI;
};

} // namespace

template <class T> Reference MakeReference(const WaveFunction::tWaveFunction<T>& wf, const OccupationConfig& occ,
                                           Reservoirs res, double eigenNoise)
{
    std::vector<ReferenceBlock> blocks;
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        ReferenceBlock b;
        b.irrep=ir;
        b.w=ir.sym->GetWeight();
        const size_t n=os.GetNumOrbitals();
        b.e.resize(n); b.f.resize(n);
        size_t i=0;
        for (const auto* o : os.Iterate())
        {
            const double g=o->GetDegeneracy();
            if (i==0) b.g=g;
            else if (g!=b.g) throw std::logic_error("Response::MakeReference: orbitals of one block with different capacities");
            b.e[i]=o->GetEigenEnergy();
            b.f[i]=o->GetOccupation()/g;
            i++;
        }
        blocks.push_back(std::move(b));
    });
    // The reservoir: one id per chemical potential, keyed on the block's k-point and spin unless the run
    // shares μ across them.
    {
        std::vector<std::string> keys;
        size_t b=0;
        for (const Irrep& ir : wf.GetQNs())
        {
            std::ostringstream key;
            if (!res.acrossK)    key << ir.sym->SequenceIndex() << ":";
            if (!res.acrossSpin) key << int(ir.ms);
            size_t id=0;
            while (id<keys.size() && keys[id]!=key.str()) id++;
            if (id==keys.size()) keys.push_back(key.str());
            blocks[b++].reservoir=int(id);
        }
    }
    return Reference(std::move(blocks), MakeOccupancyRule(occ), eigenNoise);
}
template Reference MakeReference<double>(const WaveFunction::tWaveFunction<double>&, const OccupationConfig&, Reservoirs, double);
template Reference MakeReference<dcmplx>(const WaveFunction::tWaveFunction<dcmplx>&, const OccupationConfig&, Reservoirs, double);

double EigenpairResidual(const WaveFunction::cWaveFunction& wf, Hamiltonian::cHamiltonian& H,
                         const BasisSet::tBasisSet<dcmplx>* wholeBasis)
{
    const auto D=wf.GetChargeDensity();
    double worst=0.0;
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        if constexpr (std::is_same_v<U,dcmplx>)
        {
            const auto* bs=dynamic_cast<const ChargeDensity::tobs_t<dcmplx>*>(os.GetBasisSet());
            if (!bs) throw std::logic_error("Response::EigenpairResidual: an orbital block that is not an Orbital_1E_IBS");
            const cmat_t C=Coefficients(os, bs->GetNumFunctions());
            const cmat_t CFC=blazem::ctrans(C)*cmat_t(H.GetMatrix(bs, ir.ms, D.get(), wholeBasis))*C;
            size_t i=0;
            for (const auto* o : os.Iterate())
            {
                for (size_t j=0;j<CFC.columns();j++)
                    worst=std::max(worst, std::abs(CFC(i,j)-(i==j ? dcmplx(o->GetEigenEnergy()) : dcmplx(0.0))));
                i++;
            }
        }
    });
    return worst;
}

template <class T> OrbitalFrame<T> MakeOrbitalFrame(const Reference& ref, const WaveFunction::tWaveFunction<T>& wf)
{
    std::vector<FrameBlock<T>> blocks;
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        if constexpr (!std::is_same_v<U,T>)
            throw std::logic_error("Response::MakeOrbitalFrame: a REAL TRIM block inside a complex run -- the response "
                                   "faces serve complex blocks only (the real-block sibling is not built): run the ground "
                                   "state with forceComplex");
        else
        {
            const auto* bs=dynamic_cast<const ChargeDensity::tobs_t<T>*>(os.GetBasisSet());
            if (!bs) throw std::logic_error("Response::MakeOrbitalFrame: an orbital block that is not an Orbital_1E_IBS");
            blocks.push_back({ir, bs, Coefficients(os, bs->GetNumFunctions())});
        }
    });
    return OrbitalFrame<T>(ref, std::move(blocks));
}
template OrbitalFrame<double> MakeOrbitalFrame<double>(const Reference&, const WaveFunction::tWaveFunction<double>&);
template OrbitalFrame<dcmplx> MakeOrbitalFrame<dcmplx>(const Reference&, const WaveFunction::tWaveFunction<dcmplx>&);

template <class T> OperatorProbe MakeDipoleProbe(const Reference& ref, const OrbitalFrame<T>& frame,
                                                 const WaveFunction::tWaveFunction<T>& wf, const qcMesh::Mesh& mesh)
{
    // ⚠ ONE SPATIAL BLOCK ONLY.  A dipole component that is not totally symmetric (x is B1 in C2v) couples
    // DIFFERENT irreps; under the Invariant rule its within-block matrix elements vanish by symmetry, so a
    // symmetry-adapted molecule would get a SILENTLY zero alpha_xx.  Refused until the point-group product
    // SelectionRule (and bra != ket frame pairs) exist -- doc/LinearResponsePlan.md §5c.
    for (size_t b=1;b<ref.NumBlocks();b++)
        if (ref.BlockIrrep(b).sym->SequenceIndex()!=ref.BlockIrrep(0).sym->SequenceIndex())
            throw std::logic_error("Response::MakeDipoleProbe: the reference has more than one spatial symmetry block "
                                   "(a symmetry-adapted molecule, or a k-mesh).  A dipole couples DIFFERENT irreps, which "
                                   "needs the point-group product selection rule -- run the molecule without symmetry");
    auto rule=std::make_shared<Symmetry::Invariant>();
    std::vector<Hamiltonian::AO_TransitionFock<T>> ao(3, Hamiltonian::AO_TransitionFock<T>(rule));
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        if constexpr (std::is_same_v<U,T>)
        {
            const auto* bs=dynamic_cast<const ChargeDensity::tobs_t<T>*>(os.GetBasisSet());
            if (!bs) throw std::logic_error("Response::MakeDipoleProbe: an orbital block that is not an Orbital_1E_IBS");
            for (int i=0;i<3;i++) ao[i].Add(ir, ir, mat_t<T>(qcMesh::MatrixOverlap<T>(mesh, *bs, Coordinate(i))));
        }
        else throw std::logic_error("Response::MakeDipoleProbe: a block whose scalar is not the run's");
    });
    std::vector<BlockPairs> ops;
    for (int i=0;i<3;i++) ops.push_back(frame.ToMO(ao[i], *rule));
    return OperatorProbe(ref, std::move(ops), {"x","y","z"}, rule);
}
template OperatorProbe MakeDipoleProbe<double>(const Reference&, const OrbitalFrame<double>&, const WaveFunction::tWaveFunction<double>&, const qcMesh::Mesh&);
template OperatorProbe MakeDipoleProbe<dcmplx>(const Reference&, const OrbitalFrame<dcmplx>&, const WaveFunction::tWaveFunction<dcmplx>&, const qcMesh::Mesh&);

AmplitudeProbe MakeHubbardProbe(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                                const Hamiltonian::HubbardChannels& hub)
{
    return std::move(MakeHubbardProbe(ref, wf, hub, HubbardResolution::Manifold).probe);
}

namespace {
//! An orthonormal basis of the range of the orthogonal projector \a P (its eigenvectors at eigenvalue 1).
rmat_t RangeBasis(const rmat_t& P)
{
    const size_t m=P.rows();
    rsmat_t S(m);
    for (size_t a=0;a<m;a++) for (size_t b=a;b<m;b++) S(a,b)=0.5*(P(a,b)+P(b,a));
    rvec_t lam; rmat_t v; blazem::eigen(S, lam, v);
    std::vector<size_t> keep;
    for (size_t i=0;i<m;i++)
    {
        if (std::abs(lam[i]-1.0)<1e-6) keep.push_back(i);
        else if (std::abs(lam[i])>1e-6) throw std::logic_error("Response::MakeHubbardProbe: a slot projector that is not a projector");
    }
    rmat_t Y(m, keep.size());
    for (size_t c=0;c<keep.size();c++) for (size_t a=0;a<m;a++) Y(a,c)=v(a,keep[c]);
    return Y;
}
//! "z2:1.00 x2-y2:1.00": how much of each named projector function lies in the channel spanned by \a X's columns
//! (\f$\sum_c X_{ac}^2\f$), names repeated across a multi-shell manifold averaged; weights below 0.005 left out.
//! A ONE-dimensional channel over distinct names prints its SIGNED coefficients instead ("x:+0.71 y:-0.71"), the
//! largest made positive (an eigenvector's overall sign is arbitrary): weights alone cannot tell (x+y)/sqrt2 from
//! (x-y)/sqrt2 -- on rutile O that is the sigma bond vs the plane normal (beta-MnO2's "O pz", 2026-10-22).
std::string Composition(const rmat_t& X, const std::vector<std::string>& names)
{
    if (names.size()!=X.rows()) return "";
    const bool distinct=[&]{ for (size_t a=0;a<names.size();a++) for (size_t b=a+1;b<names.size();b++) if (names[a]==names[b]) return false; return true; }();
    if (X.columns()==1 && distinct)
    {
        size_t big=0; for (size_t a=1;a<X.rows();a++) if (std::abs(X(a,0))>std::abs(X(big,0))) big=a;
        const double s = X(big,0)<0.0 ? -1.0 : 1.0;
        std::ostringstream os; os<<std::fixed<<std::setprecision(2)<<std::showpos;
        for (size_t a=0;a<X.rows();a++) if (X(a,0)*X(a,0)>=0.005) os<<(os.tellp()>0 ? " " : "")<<names[a]<<":"<<s*X(a,0);
        return os.str();
    }
    std::vector<std::string> order; std::vector<double> sum, count;
    for (size_t a=0;a<X.rows();a++)
    {
        double w=0.0; for (size_t c=0;c<X.columns();c++) w+=X(a,c)*X(a,c);
        size_t k=0; while (k<order.size() && order[k]!=names[a]) k++;
        if (k==order.size()) { order.push_back(names[a]); sum.push_back(0.0); count.push_back(0.0); }
        sum[k]+=w; count[k]+=1.0;
    }
    std::ostringstream os; os<<std::fixed<<std::setprecision(2);
    for (size_t k=0;k<order.size();k++) if (sum[k]/count[k]>=0.005) os<<(os.tellp()>0 ? " " : "")<<order[k]<<":"<<sum[k]/count[k];
    return os.str();
}
} // namespace

HubbardProbe MakeHubbardProbe(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                              const Hamiltonian::HubbardChannels& hub, HubbardResolution res)
{
    // Per block, per MANIFOLD: the amplitudes, and (for the eigen resolution) the block's weight and occupations.
    std::vector<std::vector<cmat_t>> ell;
    std::vector<double> weight;
    std::vector<rvec_t> occ;
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        const Irrep& mine=ref.BlockIrrep(ell.size());
        if (ir<mine || mine<ir)   // Irrep is ordered, not equality-comparable
            throw std::logic_error("Response::MakeHubbardProbe: the wave function's blocks are not the reference's");
        const auto* blk=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
        if (!blk) throw std::logic_error("Response::MakeHubbardProbe: an orbital block that is not an Orbital_DFT_IBS");
        std::vector<cmat_t> perManifold;
        for (const auto& l : hub.ProjectorAmplitudes(*blk, Coefficients(os, blk->GetNumFunctions())))
            perManifold.push_back(ToComplex(l));
        ell.push_back(std::move(perManifold));
        weight.push_back(ir.sym->GetWeight());
        rvec_t f(os.GetNumOrbitals()); size_t i=0;
        for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>()) f[i++]=o->GetOccupation();
        occ.push_back(std::move(f));
    });
    const auto ch=hub.Channels();
    // The channels: per manifold, the columns X (m_M x k, orthonormal) of each channel's subspace of its projector functions.
    std::vector<rmat_t> X; std::vector<size_t> man; std::vector<std::string> labels, comp;
    for (size_t M=0;M<ch.size();M++)
    {
        const size_t m=ell.empty() ? 0 : ell[0][M].rows();
        std::ostringstream base; base<<"site"<<ch[M].site<<" l="<<ch[M].l;
        const Hamiltonian::HubbardSlots hs = res==HubbardResolution::Manifold ? Hamiltonian::HubbardSlots{} : hub.SlotProjectors(M);
        if (hs.slots.empty())                                    // one channel: the manifold itself
        {
            rmat_t I(m,m,0.0); for (size_t a=0;a<m;a++) I(a,a)=1.0;
            X.push_back(I); man.push_back(M); labels.push_back(base.str()); comp.push_back("");
            continue;
        }
        rmat_t n(m,m,0.0);                                       // the spin-summed, k-summed occupation matrix (eigen only)
        if (res==HubbardResolution::Eigen)
        {
            cmat_t nc(m,m,dcmplx(0.0));
            for (size_t b=0;b<ell.size();b++)
                for (size_t i=0;i<occ[b].size();i++)
                {
                    if (occ[b][i]==0.0) continue;
                    const double f=weight[b]*occ[b][i];
                    for (size_t a=0;a<m;a++) for (size_t c=0;c<m;c++) nc(a,c)+=dcmplx(f)*ell[b][M](a,i)*std::conj(dcmplx(ell[b][M](c,i)));
                }
            for (size_t a=0;a<m;a++) for (size_t c=0;c<m;c++) n(a,c)=nc(a,c).real();
        }
        for (const auto& slot : hs.slots)
        {
            const rmat_t Y=RangeBasis(slot.P);
            if (res==HubbardResolution::Slot)
            {
                X.push_back(Y); man.push_back(M); labels.push_back(base.str()+" "+slot.label); comp.push_back(Composition(Y, hs.components));
                continue;
            }
            // Eigen: diagonalise n inside the slot, then group the eigenvalues in runs of the irrep's dimension (the partners of
            // one copy are degenerate by symmetry; the copies are not).  Ascending order, so a run is a degenerate set.
            const size_t d=Y.columns(), k=slot.irrepDim;
            if (k==0 || d%k) throw std::logic_error("Response::MakeHubbardProbe: a slot whose dimension is not a multiple of its irrep's");
            const rmat_t ns=blazem::trans(Y)*n*Y;
            rsmat_t S(d); for (size_t a=0;a<d;a++) for (size_t c=a;c<d;c++) S(a,c)=0.5*(ns(a,c)+ns(c,a));
            rvec_t lam; rmat_t Z; blazem::eigen(S, lam, Z);
            for (size_t g=0;g<d/k;g++)
            {
                rmat_t Zg(d,k);
                for (size_t c=0;c<k;c++) for (size_t a=0;a<d;a++) Zg(a,c)=Z(a,g*k+c);
                const rmat_t Xg=Y*Zg;
                double mean=0.0; for (size_t c=0;c<k;c++) mean+=lam[g*k+c]/double(k);
                std::ostringstream lab; lab<<base.str()<<" "<<slot.label<<"."<<g<<" n="<<std::fixed<<std::setprecision(3)<<mean;
                X.push_back(Xg); man.push_back(M); labels.push_back(lab.str()); comp.push_back(Composition(Xg, hs.components));
            }
        }
    }
    // Amplitudes per block per channel: X_c^T ell_M.
    std::vector<std::vector<cmat_t>> amp(ell.size());
    for (size_t b=0;b<ell.size();b++)
        for (size_t c=0;c<X.size();c++)
        {
            cmat_t Xc(X[c].rows(), X[c].columns());
            for (size_t a=0;a<X[c].rows();a++) for (size_t j=0;j<X[c].columns();j++) Xc(a,j)=X[c](a,j);
            amp[b].push_back(cmat_t(blazem::trans(Xc)*ell[b][man[c]]));
        }
    return HubbardProbe{AmplitudeProbe(ref, std::move(amp), labels), std::move(man), std::move(comp)};
}

void FeedOccupiedOrbitals(Hamiltonian::ACBN0Estimator& est, const WaveFunction::cWaveFunction& wf)
{
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        const auto* blk=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
        if (!blk) throw std::logic_error("Response::FeedOccupiedOrbitals: an orbital block that is not an Orbital_DFT_IBS");
        std::vector<const Orbitals::TOrbital<U>*> occ;
        for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>()) if (o->IsOccupied()) occ.push_back(o);
        mat_t<U> C(blk->GetNumFunctions(), occ.size());
        rvec_t f(occ.size());
        for (size_t i=0;i<occ.size();i++)
        {
            const vec_t<U>& c=occ[i]->GetCoeff();
            for (size_t a=0;a<c.size();a++) C(a,i)=c[a];
            f[i]=occ[i]->GetOccupation();
        }
        est.Accumulate(*blk, ir.ms, ir.sym->GetWeight(), C, f);
    });
}

} // namespace

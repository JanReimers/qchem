// File: Response/Imp/Response.C  Walk a converged Bloch wave function into the response vocabulary.
module;
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Response;
import qchem.Orbitals;                   // TOrbitals / TOrbital (the per-block orbital sets)
import qchem.BasisSet.Orbital_DFT_IBS;   // the block basis the projector amplitudes are asked on
import qchem.Blaze;

namespace qchem::Response
{

namespace {

// Every block is visited in the wave function's own QN order, by both adapters -- so block b of the
// Reference and block b of the probe are the same (k, σ).  A block's orbitals are either scalar (a real
// TRIM block inside a complex run is TOrbitals<double>); the visitor hands each branch its native type.
template <class Visit> void ForEachBlock(const WaveFunction::cWaveFunction& wf, Visit&& visit)
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

} // namespace

Reference MakeReference(const WaveFunction::cWaveFunction& wf, const OccupationConfig& occ, Reservoirs res,
                        double eigenNoise)
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

AmplitudeProbe MakeHubbardProbe(const Reference& ref, const WaveFunction::cWaveFunction& wf,
                                const Hamiltonian::HubbardChannels& hub)
{
    std::vector<std::vector<cmat_t>> amp;
    ForEachBlock(wf, [&]<class U>(const Irrep& ir, const Orbitals::TOrbitals<U>& os)
    {
        const Irrep& mine=ref.BlockIrrep(amp.size());
        if (ir<mine || mine<ir)   // Irrep is ordered, not equality-comparable
            throw std::logic_error("Response::MakeHubbardProbe: the wave function's blocks are not the reference's");
        const auto* blk=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(os.GetBasisSet());
        if (!blk) throw std::logic_error("Response::MakeHubbardProbe: an orbital block that is not an Orbital_DFT_IBS");
        std::vector<cmat_t> perChannel;
        for (const auto& l : hub.ProjectorAmplitudes(*blk, Coefficients(os, blk->GetNumFunctions())))
            perChannel.push_back(ToComplex(l));
        amp.push_back(std::move(perChannel));
    });
    std::vector<std::string> labels;
    for (const auto& c : hub.Channels())
    {
        std::ostringstream os;
        os << "site" << c.site << " l=" << c.l;
        labels.push_back(os.str());
    }
    return AmplitudeProbe(ref, std::move(amp), std::move(labels));
}

} // namespace

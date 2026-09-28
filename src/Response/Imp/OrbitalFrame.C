// File: Response/Imp/OrbitalFrame.C  The AO <-> MO bridge.
module;
#include <cmath>
#include <complex>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>
module qchem.Response.OrbitalFrame;
import qchem.Streamable;   // Irrep's operator<<
import qchem.Blaze;

namespace qchem::Response
{

namespace {

template <class T> cmat_t Complexify(const mat_t<T>& m)
{
    cmat_t c(m.rows(), m.columns());
    for (size_t i=0;i<m.rows();i++) for (size_t j=0;j<m.columns();j++) c(i,j)=dcmplx(m(i,j));
    return c;
}

//! The Hermitian part of \a M as a \c hmat_t<T>, THROWING if \a M is not Hermitian (relative 1e-8), or -- on a
//! real block -- not real.  Neither is a rounding question: R0 maps Hermitian to Hermitian, so either is a defect.
template <class T> hmat_t<T> HermitianAO(const cmat_t& M, const Irrep& ir)
{
    double scale=0, antiH=0, imag=0;
    for (size_t i=0;i<M.rows();i++)
        for (size_t j=0;j<M.columns();j++)
        {
            scale=std::max(scale, std::abs(M(i,j)));
            antiH=std::max(antiH, std::abs(M(i,j)-std::conj(M(j,i))));
            imag =std::max(imag,  std::fabs(M(i,j).imag()));
        }
    const double tol=1e-8*scale+1e-300;
    std::ostringstream os;
    if (antiH>tol) os << "OrbitalFrame::ToAO: block " << ir << " -- δD is not Hermitian (" << antiH << " vs " << scale << ")";
    if constexpr (std::is_floating_point_v<T>)
        if (imag>tol) os << "OrbitalFrame::ToAO: block " << ir << " is REAL but δD has an imaginary part (" << imag << ")";
    if (!os.str().empty()) throw std::logic_error(os.str());
    const size_t n=M.rows();
    hmat_t<T> H(n);
    for (size_t i=0;i<n;i++)
        for (size_t j=i;j<n;j++)
        {
            const dcmplx h=0.5*(M(i,j)+std::conj(M(j,i)));
            if constexpr (std::is_floating_point_v<T>) H(i,j)=h.real();
            else                                       H(i,j)=h;
        }
    return H;
}

} // namespace

template <class T> OrbitalFrame<T>::OrbitalFrame(const Reference& ref, std::vector<FrameBlock<T>> blocks)
    : itsRef(ref), itsBlocks(std::move(blocks))
{
    if (itsBlocks.size()!=itsRef.NumBlocks()) throw std::invalid_argument("OrbitalFrame: not one frame block per reference block");
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        const auto& f=itsBlocks[b];
        const Irrep& r=itsRef.BlockIrrep(b);
        if (f.irrep<r || r<f.irrep) throw std::invalid_argument("OrbitalFrame: the blocks are not in the reference's order");
        if (!f.bs) throw std::invalid_argument("OrbitalFrame: a block with no basis");
        if (f.C.columns()!=itsRef.NumOrbitals(b) || f.C.rows()!=f.bs->GetNumFunctions())
            throw std::invalid_argument("OrbitalFrame: a block's coefficients do not match its basis / orbital count");
    }
}

template <class T> std::unique_ptr<ChargeDensity::TransitionDensity<T>>
OrbitalFrame<T>::ToAO(const BlockPairs& dD, std::shared_ptr<const Symmetry::SelectionRule> rule) const
{
    const std::vector<size_t> p=itsRef.Partners(*rule);
    if (dD.m.size()!=itsBlocks.size()) throw std::invalid_argument("OrbitalFrame::ToAO: one matrix per block");
    std::vector<ChargeDensity::TransitionBlock<T>> out;
    out.reserve(itsBlocks.size());
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        if (p[b]!=b) throw std::logic_error("OrbitalFrame::ToAO: a (bra != ket) block pair -- the AO transition density "
                                            "of a q != 0 / symmetry-lowering perturbation is stage R3's");
        const cmat_t C=Complexify(itsBlocks[b].C);
        const cmat_t M=C*dD.m[b]*blazem::ctrans(C);
        out.push_back({itsBlocks[b].irrep, itsBlocks[b].bs, HermitianAO<T>(M, itsBlocks[b].irrep)});
    }
    return ChargeDensity::AO_TransitionDensity_Factory<T>(std::move(out), std::move(rule));
}

template <class T> BlockPairs OrbitalFrame<T>::ToMO(const Hamiltonian::TransitionFock<T>& dF, const Symmetry::SelectionRule& rule) const
{
    const std::vector<size_t> p=itsRef.Partners(rule);
    BlockPairs X;
    X.m.reserve(itsBlocks.size());
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        const auto& ket=itsBlocks[b];
        const auto& bra=itsBlocks[p[b]];
        const cmat_t F=Complexify(dF.Matrix(bra.irrep, ket.irrep));
        X.m.push_back(cmat_t(blazem::ctrans(Complexify(bra.C))*F*Complexify(ket.C)));
    }
    return X;
}

template class OrbitalFrame<double>;
template class OrbitalFrame<dcmplx>;

} // namespace

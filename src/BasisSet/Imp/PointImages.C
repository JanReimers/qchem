// File: BasisSet/Imp/PointImages.C  The image record's bodies: building it, and forming a batch of Bloch rows from it.
module;
#include <cstddef>
#include <complex>     // dcmplx * double (the phase sum)
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
module qchem.BasisSet.PointImages;

namespace qchem::BasisSet
{

PointImages::PointImages(size_t nInner, std::vector<ivec3_t> offsets)
    : itsNInner(nInner), itsNOut(nInner), itsOffsets(std::move(offsets))
{
    if (itsOffsets.size()>size_t(UINT32_MAX) || nInner>size_t(UINT32_MAX))
        throw std::length_error("PointImages: more images or functions than a 32-bit index holds");
}

void PointImages::AddImage(size_t img, const rvec_t& chi)
{
    if (img>=itsOffsets.size()) throw std::out_of_range("PointImages::AddImage: image index past the offset list");
    if (chi.size()!=itsNInner)  throw std::invalid_argument("PointImages::AddImage: an image row is not one value per inner function");
    bool any=false;
    for (size_t c=0; c<chi.size(); c++)
        if (chi[c]!=0.0)
        {
            itsCol.push_back(uint32_t(c));
            itsVal.push_back(chi[c]);
            any=true;
        }
    if (!any) return;                      // an all-zero image adds nothing (and +0 accumulators stay put)
    itsSegImg.push_back(uint32_t(img));
    itsSegBegin.push_back(itsVal.size());
}

void PointImages::EndPoint() {itsRow.push_back(itsSegImg.size());}

void PointImages::Append(const PointImages& t)
{
    if (t.itsNInner!=itsNInner || t.itsOffsets.size()!=itsOffsets.size())
        throw std::invalid_argument("PointImages::Append: the chunks are over different functions or offsets");
    if (t.itsNOut!=itsNOut || t.itsTransform!=itsTransform)   // a transformed view's chunks each carry the view's one transform
        throw std::invalid_argument("PointImages::Append: the chunks carry different transforms");
    const size_t seg0=itsSegImg.size(), ent0=itsVal.size();
    for (size_t g=1; g<t.itsRow.size(); g++)      itsRow.push_back(seg0+t.itsRow[g]);
    for (size_t j=0; j<t.itsSegImg.size(); j++)
    {
        itsSegImg.push_back(t.itsSegImg[j]);
        itsSegBegin.push_back(ent0+t.itsSegBegin[j+1]);
    }
    itsCol.insert(itsCol.end(), t.itsCol.begin(), t.itsCol.end());
    itsVal.insert(itsVal.end(), t.itsVal.begin(), t.itsVal.end());
}

void PointImages::SetTransform(Transform t, size_t nOut)
{
    if (t.size()!=nOut) throw std::invalid_argument("PointImages::SetTransform: one coefficient list per output function");
    for (const auto& col : t)
        for (const auto& [c,v] : col)
            if (c>=itsNInner) throw std::out_of_range("PointImages::SetTransform: an inner function index past the inner count");
    itsTransform=std::move(t);
    itsNOut=nOut;
}

double PointImages::Bytes() const
{
    double b=double(itsVal.size())*sizeof(double) + double(itsCol.size())*sizeof(uint32_t)
            +double(itsSegImg.size())*sizeof(uint32_t) + double(itsSegBegin.size()+itsRow.size())*sizeof(size_t)
            +double(itsOffsets.size())*sizeof(ivec3_t);
    for (const auto& col : itsTransform) b+=double(col.size())*sizeof(std::pair<size_t,double>);
    return b;
}

// The dense build's order exactly: per point, images ascending, Phi(q,c) += phase*chi (dcmplx * double), then the transform per point.
mat_t<dcmplx> PointImages::Bloch(const cvec_t& phases, size_t lo, size_t hi) const
{
    if (phases.size()!=itsOffsets.size()) throw std::invalid_argument("PointImages::Bloch: one phase per offset, got "+std::to_string(phases.size())
                                                                       +" for "+std::to_string(itsOffsets.size()));
    if (lo>hi || hi>NumPoints()) throw std::out_of_range("PointImages::Bloch: row range past the point set");
    const size_t nb=hi-lo;
    mat_t<dcmplx> C(nb, itsNInner, dcmplx(0.0));
    for (size_t q=0; q<nb; q++)
        for (size_t j=itsRow[lo+q]; j<itsRow[lo+q+1]; j++)
        {
            const dcmplx ph=phases[itsSegImg[j]];
            for (size_t e=itsSegBegin[j]; e<itsSegBegin[j+1]; e++) C(q,itsCol[e])+=ph*itsVal[e];
        }
    if (itsTransform.empty()) return C;
    mat_t<dcmplx> P(nb, itsNOut, dcmplx(0.0));
    for (size_t q=0; q<nb; q++)
        for (size_t s=0; s<itsNOut; s++)
        {
            dcmplx a(0.0);
            for (const auto& [c,t] : itsTransform[s]) a+=t*C(q,c);
            P(q,s)=a;
        }
    return P;
}

} // namespace qchem::BasisSet

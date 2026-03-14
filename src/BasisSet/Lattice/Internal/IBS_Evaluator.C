// File: BasisSet/Lattice/Internal/IEClient.C Store basis set data needed for integral engines.
module;
#include <valarray>
#include <complex>
export module qchem.BasisSet.Lattice.Internal.IBS_Evaluator;
export import qchem.VectorFunction;
export import qchem.Types;
import qchem.Conversions;
import oml.Vector;

std::valarray<RVec3> torvec(const std::valarray<IVec3> i3)
{
    std::valarray<RVec3> r3(i3.size());
    size_t n=0;
    for (auto i:i3) r3[n++]=i;
    return r3;
}


export namespace PlaneWave
{

struct IBS_Evaluator : public virtual VectorFunction<dcmplx>
{
    using Vec    =VectorFunction<dcmplx>::Vec;
    using Vec3Vec=VectorFunction<dcmplx>::Vec3Vec;

    IBS_Evaluator(RVec3 _k,const std::valarray<IVec3>& _Gs,double _norm)
        : k(_k), Gs(torvec(_Gs)), norm(_norm)
        {}

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual size_t size() const {return Gs.size();}


    RVec3 k;
    std::valarray<RVec3> Gs;
    double norm;

};

IBS_Evaluator::Vec IBS_Evaluator::operator() (const RVec3& r) const
{
    std::valarray<dcmplx> ephi(size());
    size_t n=0;
    for (auto G:Gs) ephi[n++]=exp(dcmplx(0.0,(k+G)*r));

    return to_omlVector(ephi);
}

} //namespace
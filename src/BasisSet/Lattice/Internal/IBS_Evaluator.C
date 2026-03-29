// File: BasisSet/Lattice/Internal/IEClient.C Store basis set data needed for integral engines.
module;
#include <valarray>
#include <complex>
export module qchem.BasisSet.Lattice.Internal.IBS_Evaluator;
export import qchem.VectorFunction;
export import qchem.Types;
import qchem.Conversions;

// std::valarray<rvec3_t> torvec(const std::valarray<ivec3_t> i3)
// {
//     std::valarray<rvec3_t> r3(i3.size());
//     size_t n=0;
//     for (auto i:i3) r3[n++]=i;
//     return r3;
// }


export namespace PlaneWave
{

struct IBS_Evaluator : public virtual VectorFunction<dcmplx>
{

    IBS_Evaluator(rvec3_t _k,const std::valarray<ivec3_t>& _Gs,double _norm)
        : k(_k), Gs(torvec(_Gs)), norm(_norm)
        {}

    virtual vec_t    <dcmplx> operator() (const rvec3_t&) const;
    virtual vec3vec_t<dcmplx> Gradient   (const rvec3_t&) const;
    virtual size_t size() const {return Gs.size();}


    rvec3_t k;
    std::valarray<rvec3_t> Gs;
    double norm;

};

IBS_Evaluator::Vec IBS_Evaluator::operator() (const rvec3_t& r) const
{
    std::valarray<dcmplx> ephi(size());
    size_t n=0;
    for (auto G:Gs) ephi[n++]=exp(dcmplx(0.0,(k+G)*r));

    return to_omlVector(ephi);
}

} //namespace
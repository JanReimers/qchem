// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <valarray>
#include <vector>
export module BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Atom.Internal.ExponentGrouper;
export import qchem.BasisSet.Atom.internal.Rk;
export import qchem.BasisSet.Internal.ERI3;
export import qchem.Fit_IBS;
export import qchem.VectorFunction;
export import oml;

export using dERI3=ERI3<double>;

export template <class T> Vector<T> convert(const std::valarray<T>& v) 
{
    Vector<T> ret(v.size());
    size_t i=0;
    for (auto vi:v) ret(++i)=vi;
    return ret;
}
export template <class T> std::valarray<T> convert(const Vector<T>& v) 
{
    std::valarray<T> ret(v.size());
    size_t i=0;
    for (auto vi:v) ret[i++]=vi;
    return ret;
}

export class IBS_Evaluator : public VectorFunction<double>
{
public:
    using ds_t=std::valarray<double>;
    using is_t=std::vector<int>;
    using omls_t=SMatrix<double>;
    using omlm_t= Matrix<double>;
    using omlv_t= Vector<double>;
    virtual ~IBS_Evaluator() {};

    virtual void Register(ExponentGrouper*)=0; //Set up unique spline or exponent indexes.
    virtual size_t size() const =0;
    virtual int    Getl() const =0;

    virtual omls_t Overlap   () const=0;
    virtual omls_t Grad2     () const=0;
    virtual omls_t Inv_r1    () const=0;
    virtual omls_t Inv_r2    () const=0;
    virtual omls_t Repulsion () const=0;
    virtual omlv_t Charge    () const=0;
    virtual ds_t   Norm      () const=0;
    virtual omlm_t XRepulsion(const Fit_IBS&) const=0;

    virtual dERI3  Overlap  (const IBS_Evaluator*) const=0; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator*) const=0; //3 center
    virtual Rk*    CreateRk (size_t ia,size_t ic,size_t ib,size_t id) const=0; //4 center


};
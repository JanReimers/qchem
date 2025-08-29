// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <valarray>
#include <vector>
#include <ranges>
#include <iosfwd>
export module BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Atom.Internal.ExponentGrouper;
export import qchem.BasisSet.Atom.internal.Rk;
export import qchem.BasisSet.Internal.ERI3;
export import qchem.Fit_IBS;
export import qchem.Orbital_DHF_IBS;
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
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
public:
    using ds_t=std::valarray<double>;
    using is_t=std::vector<int>;
    using omls_t=SMatrix<double>;
    using omlm_t= Matrix<double>;
    using omlv_t= Vector<double>;

    IBS_Evaluator(int _l, const is_t& _mls) : l(_l), mls(_mls),ns(0),grouper(0) {};
    virtual ~IBS_Evaluator() {};

    virtual void          Register(Grouper*)=0; //Set up unique spline or exponent indexes.
    virtual size_t        size    (             ) const {return ns.size();}
    virtual int           Getl    (             ) const {return l;}
    virtual size_t        es_index(size_t i     ) const {return es_indices[i];}
    virtual const is_t&   Getmls  (             ) const {return mls;}

    iota_view             indices (             ) const {return iota_view(size_t(0),size());}
    iota_view             indices (size_t start ) const {return iota_view(start,size());}
    virtual std::ostream& Write   (std::ostream&) const=0;
    virtual size_t maxSpan() const {return size();}  //assume no overlap for indeces separated by > maxSpan

    virtual size_t        GetVectorSize() const {return size();}

    virtual omls_t Overlap   () const=0;
    virtual omls_t Grad2     () const=0;
    virtual omls_t Inv_r1    () const=0;
    virtual omls_t Inv_r2    () const=0;
    virtual omls_t Repulsion () const=0;
    virtual omlv_t Charge    () const=0;
    virtual ds_t   Norm      () const=0;
    virtual omlm_t XRepulsion(const Fit_IBS&) const=0;
    virtual omlm_t XKinetic  (const Orbital_RKBS_IBS<double>*) const=0;

    virtual dERI3  Overlap  (const Fit_IBS&) const=0; //3 center
    virtual dERI3  Repulsion(const Fit_IBS&) const=0; //3 center
protected:
    int  l;
    is_t mls;
    ds_t ns;
    const ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index

};
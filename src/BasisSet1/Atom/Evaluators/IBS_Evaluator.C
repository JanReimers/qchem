// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <vector>
#include <ranges>
#include <iosfwd>
#include <cassert>
#include <sstream>
export module qchem.BasisSet1.Atom.Evaluators.IBS;
export import qchem.BasisSet1.Atom.Evaluators.Internal.ExponentGrouper;
export import qchem.BasisSet1.Internal.ERI3;
export import qchem.Symmetry.Irrep;
export import qchem.VectorFunction;
import qchem.Symmetry.Ylm;

export using dERI3=ERI3<double>;

export template <class E> concept isEvaluator = requires (E e,size_t i, size_t j, size_t ic)
            {
                e.Overlap(i,j); //Should all be inline.
                e.Grad2  (i,j);
                e.Inv_r1 (i,j);
                e.Inv_r2 (i,j);
                e.Overlap  (i,j,e,ic); 
                e.Repulsion(i,j,e,ic);
                e.Repulsion(i,j);
                e.Charge(i);
                e.Norm(i);
                // etc.....
            };

export class IBS_Evaluator : public VectorFunction<double>
{
    typedef std::ranges::iota_view<size_t,size_t> iota_view;
public:
    using is_t=std::vector<int>;
    

    IBS_Evaluator(int _l, const is_t& _mls) : l(_l), mls(_mls),ns(0),grouper(0) {};
    IBS_Evaluator(const Irrep_QNs::sym_t& ylm);
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

    virtual rsmat_t Repulsion () const=0;
    virtual  rvec_t Charge    () const=0;
    virtual  rvec_t Norm      () const=0;
    
    virtual rmat_t XRepulsion(const IBS_Evaluator&) const=0;
    virtual rmat_t XKinetic  (const IBS_Evaluator&) const=0;
    virtual std::string RadialID () const=0;
    virtual std::string AngularID() const;
    virtual std::string Name() const=0;
protected:
    int  l;
private:
    // used by Coulomb_AngularIntegrals and ExchangeAngularIntegrals
    is_t   mls;
protected:
    rvec_t ns;
    const ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index

};



IBS_Evaluator::IBS_Evaluator(const Irrep_QNs::sym_t& y) :  l(0), mls({}),ns(0),grouper(0)
{
    const Yl_Sym* yl=dynamic_cast<const Yl_Sym*>(y.get());
    assert(yl);
    l=yl->GetL();
    const Ylm_Sym* ylm=dynamic_cast<const Ylm_Sym*>(y.get());
    if (ylm)
    {
        mls=ylm->Getmls();
    }
}

std::string IBS_Evaluator::AngularID() const
{
     std::ostringstream os;
     os << l << " {";
     for (auto ml:mls) os << ml << " ";
     os << "}";
     return os.str();
}
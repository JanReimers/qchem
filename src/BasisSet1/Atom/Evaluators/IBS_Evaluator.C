// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <vector>
#include <ranges>
#include <iosfwd>
#include <cassert>
#include <sstream>
export module qchem.BasisSet1.Atom.Evaluators.IBS;
export import qchem.BasisSet1.Atom.Internal.ExponentGrouper;
export import qchem.BasisSet1.Internal.ERI3;
// export import qchem.Fit_IBS;
// export import qchem.Orbital_DHF_IBS;
export import qchem.Symmetry.Irrep;
export import qchem.VectorFunction;
import qchem.Symmetry.Ylm;

export using dERI3=ERI3<double>;

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

    virtual rsmat_t Overlap   () const=0;
    virtual rsmat_t Grad2     () const=0;
    virtual rsmat_t Inv_r1    () const=0;
    virtual rsmat_t Inv_r2    () const=0;
    virtual rsmat_t Repulsion () const=0;
    virtual  rvec_t Charge    () const=0;
    virtual  rvec_t Norm      () const=0;
    virtual rmat_t XRepulsion(const IBS_Evaluator&) const=0;
    virtual rmat_t XKinetic(const IBS_Evaluator*) const=0;

    // virtual rmat_t XRepulsion(const Fit_IBS& b) const
    // {
    //     return XRepulsion(dynamic_cast<const IBS_Evaluator&>(b));
    // }

    // virtual rmat_t XKinetic  (const Orbital_RKBS_IBS<double>*b) const
    // {
    //     return XKinetic(dynamic_cast<const IBS_Evaluator*>(b));
    // }

    // virtual dERI3  Overlap  (const Fit_IBS& c) const //3 center
    // {
    //     return Overlap(dynamic_cast<const IBS_Evaluator&>(c));
    // }
    // virtual dERI3  Repulsion(const Fit_IBS& c) const //3 center
    // {
    //     return Repulsion(dynamic_cast<const IBS_Evaluator&>(c));
    // }
    virtual dERI3  Overlap  (const IBS_Evaluator&) const=0; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator&) const=0; //3 center
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
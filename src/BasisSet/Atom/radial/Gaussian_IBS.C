// File: BasisSet/Atom/radial/Gaussian_IBS.C
module;
#include <valarray>
#include <vector>
export module BasisSet.Atom.Gaussian_IBS;
import BasisSet.Atom.IBS_Evaluator;

export class Gaussian_IBS : public virtual IBS_Evaluator
{
public: 
 
    Gaussian_IBS(const ds_t& _es, int _l, const is_t& _mls) : es(_es), l(_l), mls(_mls),ns(norms()), grouper(0) {};
    Gaussian_IBS(const omlv_t& _es, int _l, const is_t& _mls) : es(convert(_es)), l(_l), mls(_mls),ns(norms()), grouper(0) {};
    virtual void Register(ExponentGrouper*); //Set up unique spline or exponent indexes.
    virtual size_t size() const {return es.size();}
    virtual int    Getl() const {return l;}
    virtual size_t GetVectorSize() const {return size();}

    virtual omls_t Overlap  () const;
    virtual omls_t Grad2    () const;
    virtual omls_t Inv_r1   () const;
    virtual omls_t Inv_r2   () const;
    virtual omls_t Repulsion() const;
    virtual omlv_t Charge   () const;
    virtual ds_t   Norm     () const {return ns;}

    virtual dERI3  Overlap  (const IBS_Evaluator*) const; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator*) const; //3 center
    
    virtual Rk* CreateRk(size_t ia,size_t ic,size_t ib,size_t id) const;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

protected:
    ds_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;

    ds_t es; 
    int  l;
    is_t mls;
    ds_t ns;
    const ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index
};


export class Gaussian_RKBS_IBS : public Gaussian_IBS
{
    public:
    Gaussian_RKBS_IBS(const omlv_t& es, int l, const is_t& mls) : Gaussian_IBS(es,l,mls) {ns=norms();}
    virtual ds_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;
};
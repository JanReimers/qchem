// File: BasisSet/Atom/radial/Slater_IBS.C
module;
#include <valarray>
#include <vector>
#include <iosfwd>
export module BasisSet.Atom.Slater_IBS;
export import BasisSet.Atom.Exponential_IBS_Evaluator;

export class Slater_IBS : public Exponential_IBS_Evaluator
{
public: 
 
    Slater_IBS(const   ds_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}
    Slater_IBS(const omlv_t& es, int l, const is_t& mls) : Exponential_IBS_Evaluator(es,l,mls) {ns=norms();}

   
    virtual std::ostream& Write   (std::ostream&) const;
 
    virtual omls_t Overlap  () const;
    virtual omls_t Grad2    () const;
    virtual omls_t Inv_r1   () const;
    virtual omls_t Inv_r2   () const;
    virtual omls_t Repulsion() const;
    virtual omlv_t Charge   () const;
    virtual ds_t   Norm     () const {return ns;}
    virtual omlm_t XRepulsion(const Fit_IBS&) const;
    virtual omlm_t XKinetic  (const Orbital_RKBS_IBS<double>*) const;

    virtual dERI3  Overlap  (const Fit_IBS&) const; //3 center
    virtual dERI3  Repulsion(const Fit_IBS&) const; //3 center

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

protected:
    ds_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;

   
};

export class Slater_RKBS_IBS : public Slater_IBS
{
public:
    Slater_RKBS_IBS(const omlv_t& es, int _kappa, int l,const is_t& mls) : Slater_IBS(es,l,mls), kappa(_kappa) {ns=norms();};
    ds_t norms() const; //assumes es,l are already initialized
    virtual double Inv_r1(double ea , double eb,size_t l_total) const;

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
private:
    ds_t eval(const RVec3&) const;
    int kappa;
};
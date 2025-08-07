// File: BasisSet/Atom/radial/Slater_IBS.C
module;
#include <valarray>
#include <vector>
export module BasisSet.Atom.Slater_IBS;
import BasisSet.Atom.IBS_Evaluator;

export class Slater_IBS : public virtual IBS_Evaluator
{
public: 
 
    Slater_IBS(const ds_t& _es, int _l, const is_t& _mls) : es(_es), l(_l), mls(_mls),ns(norms()) {};
    virtual void Register(ExponentGrouper&); //Set up unique spline or exponent indexes.

    virtual size_t size() const {return es.size();}

    virtual omls_t Overlap  () const;
    virtual omls_t Grad2    () const;
    virtual omls_t Inv_r1   () const;
    virtual omls_t Inv_r2   () const;
    virtual omls_t Repulsion() const;
    virtual omlv_t Charge   () const;

    // virtual ERI3   Overlap  (const Slater_IBS&) const; //3 center
    // virtual ERI3   Repulsion(const Slater_IBS&) const; //3 center

    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;

private:
    ds_t norms() const; //assumes es,l are already initialized

    ds_t es; 
    int  l;
    is_t mls;
    ds_t ns;
    std::vector<size_t> es_indices; //Unique exponent index
};
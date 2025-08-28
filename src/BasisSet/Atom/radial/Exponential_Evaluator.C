// File: BasisSet/Atom/radial/Exponential_Evaluator.C  COmmon for Slater and Gaussian evaluators
module;
#include <valarray>
#include <vector>
#include <iosfwd>

export module BasisSet.Atom.Exponential_IBS_Evaluator;

export import BasisSet.Atom.IBS_Evaluator;

export class Exponential_IBS_Evaluator : public virtual IBS_Evaluator
{
public:
    Exponential_IBS_Evaluator(const   ds_t& _es, int _l, const is_t& _mls) : es(        _es ), l(_l), mls(_mls),ns(es.size()),grouper(0) {};
    Exponential_IBS_Evaluator(const omlv_t& _es, int _l, const is_t& _mls) : Exponential_IBS_Evaluator(convert(_es),_l,_mls) {};

    virtual size_t        size    (             ) const {return es.size();}
    virtual int           Getl    (             ) const {return l;}
    virtual size_t        es_index(size_t i     ) const {return es_indices[i];}
    virtual const is_t&   Getmls  (             ) const {return mls;}
protected:
    ds_t es; 
    int  l;
    is_t mls;
    ds_t ns;
    const ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index
};
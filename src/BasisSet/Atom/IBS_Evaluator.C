// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <valarray>
#include <vector>
export module BasisSet.Atom.IBS_Evaluator;
export import qchem.VectorFunction;
export import oml;


export class IBS_Evaluator : public VectorFunction<double>
{
public:
    using ds_t=std::valarray<double>;
    using is_t=std::vector<int>;
    using omls_t=SMatrix<double>;
    using omlm_t= Matrix<double>;
    using omlv_t= Vector<double>;

    virtual size_t size() const =0;

    virtual omls_t Overlap  () const=0;
    virtual omls_t Grad2    () const=0;
    virtual omls_t Inv_r1   () const=0;
    virtual omls_t Inv_r2   () const=0;
    virtual omls_t Repulsion() const=0;
    virtual omlv_t Charge   () const=0;

};
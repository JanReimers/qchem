// File: BasisSet/Atom/IBS_Evaluator.C
module;
#include <valarray>
#include <vector>
export module BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Atom.Internal.ExponentGrouper;
export import qchem.BasisSet.Atom.internal.Rk;
export import qchem.VectorFunction;
export import oml;

export template <class T> using ERI3=std::vector<SMatrix<T>>;
export using dERI3=ERI3<double>;

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

    virtual omls_t Overlap  () const=0;
    virtual omls_t Grad2    () const=0;
    virtual omls_t Inv_r1   () const=0;
    virtual omls_t Inv_r2   () const=0;
    virtual omls_t Repulsion() const=0;
    virtual omlv_t Charge   () const=0;

    virtual dERI3  Overlap  (const IBS_Evaluator*) const=0; //3 center
    virtual dERI3  Repulsion(const IBS_Evaluator*) const=0; //3 center
    virtual Rk*    CreateRk (size_t ia,size_t ic,size_t ib,size_t id) const=0; //4 center


};
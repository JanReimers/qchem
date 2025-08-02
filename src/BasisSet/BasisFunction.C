// File: BasisFunction.H
export module qchem.BasisFunction;
export import qchem.ScalarFunction;
import qchem.Streamable;

//---------------------------------------------------------------
//
//  Basis function is nothing more than a scalar funtion op(r) 
//  that can be output streamed and cloned.
//
export template <class T> class TBasisFunction
    : public virtual Streamable
    , public virtual ScalarFunction<T>
{
public:
    virtual ~TBasisFunction()  {};
    virtual TBasisFunction* Clone() const=0;
};

export using    Real_BF=TBasisFunction<double>;
export using Complex_BF=TBasisFunction<dcmplx>;

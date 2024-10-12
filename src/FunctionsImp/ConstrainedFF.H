// File: ConstrainedFF.H  General constrained fit.
#ifndef _ConstrainedFF_H_
#define _ConstrainedFF_H_

#include "FunctionsImp/FittedFunctionImplementation.H"

template <class T> class ConstrainedFF
    : public FittedFunctionImplementation<T>
{
    typedef FittedFunctionImplementation<T> Base;
    typedef typename Base::Vec Vec;
public:
    ConstrainedFF();
    ConstrainedFF(const rc_ptr<IrrepBasisSet>&, const Vec&, Mesh* m);

    virtual double DoFit(const ScalarFFClient&);
    virtual double DoFit(const DensityFFClient&);

    virtual std::ostream& Write    (std::ostream&) const;
    virtual std::istream& Read     (std::istream&)      ;

private:
    using FittedFunctionImplementation<T>::itsInvRepl;
    Vec g,gS;
    T   gSg;
};

#endif //_ConstrainedFF_H_


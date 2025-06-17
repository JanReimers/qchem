// File: ConstrainedFF.H  constrained fit.
#ifndef _ConstrainedFF_H_
#define _ConstrainedFF_H_

#include "FittedFunction.H"

template <class T> class ConstrainedFF
    : public FittedFunctionImp<T>
{
    typedef FittedFunctionImp<T> Base;
    typedef typename Base::Vec Vec;
public:
    typedef typename Base::mesh_t mesh_t;
    typedef typename Base::bs_t   bs_t;

    ConstrainedFF();
    ConstrainedFF(bs_t&, const Vec&, mesh_t&  m);

    virtual double DoFit(const ScalarFFClient&);
    virtual double DoFit(const DensityFFClient&);

    virtual std::ostream& Write    (std::ostream&) const;
    virtual std::istream& Read     (std::istream&)      ;

private:
    using Base::itsLAParams;
    Vec g,gS;
    T   gSg;
};

#endif //_ConstrainedFF_H_


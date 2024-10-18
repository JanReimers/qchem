// File: SlaterExchange.H Slater exchange potential.
#ifndef _SlaterExchange_H_
#define _SlaterExchange_H_



#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Misc/Spin.H"

class SlaterExchange
    : public  virtual ExchangeFunctional
    , private ExchangeFunctionalImplementation
{
    typedef Vector3D<double> Vec3;
public:
    SlaterExchange(               );
    SlaterExchange(double theAlpha);
    SlaterExchange(double theAlpha, const Spin&);

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;
    virtual double GetVxc(double ChargeDensity) const;


    virtual std::ostream&            Write(std::ostream&) const;
    virtual std::istream&            Read (std::istream&)      ;

private:
    double itsAlpha;
    Spin   itsSpin;
};

#endif //_SlaterExchange_H_

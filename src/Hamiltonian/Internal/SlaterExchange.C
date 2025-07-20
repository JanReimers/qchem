// File: SlaterExchange.C Slater exchange potential.
module;
#include <iosfwd>
export module qchem.Hamiltonian.Internal.SlaterExchange;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Symmetry.Spin;
import oml.Vector3D;

export class SlaterExchange
    : public  ExFunctional
{
    typedef Vector3D<double> Vec3;
public:
    SlaterExchange(               );
    SlaterExchange(double theAlpha);
    SlaterExchange(double theAlpha, const Spin&);

    virtual double operator()(const Vec3&) const;
    virtual Vec3   Gradient  (const Vec3&) const;
    virtual double GetVxc(double ChargeDensity) const;


    virtual std::ostream& Write(std::ostream&) const;

private:
    double itsAlpha;
    Spin   itsSpin;
};


// File: ExchangeFunctional.C   Exchange potential for DFT.
export module qchem.Hamiltonian.Internal.ExFunctional;
import qchem.ChargeDensity;
import qchem.ScalarFunction;
import qchem.Streamable;

export class ExFunctional
    : public virtual Streamable
    , public virtual ScalarFunction<double>{
public:
    ExFunctional(               );

    virtual void           InsertChargeDensity(const DM_CD*);
    virtual Vector<double> GetVxcs(const Vector<double>& ChargeDensities) const;
    virtual double         GetVxc(                double ChargeDensity) const=0;
    virtual void           SetPolarized(bool p) {isPolarized=p;}

protected:

    const DM_CD* itsChargeDensity;
    bool            isPolarized;
};


// File: SCFIterator.H  Interface for an object that manages SCF convergence.



class ChargeDensity;
class FittedCD;
class WaveFunction;

class SCFIterator
{
public:
    virtual ~SCFIterator() {};

    virtual bool   Iterate        (double relax, double epsRo, int Nmax, double Smear)=0;
    virtual void   DumpElectrons  (WaveFunction*, double Smear)=0;
    virtual double DisplayEnergies(int i, double lam, double chargeDensityChange, double fitError) const=0;
    virtual void   DisplayEigen   () const=0;

    virtual ChargeDensity* GetExactChargeDensity () const=0;
};


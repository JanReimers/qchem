// File:  AuxillaryFJ.C  Class for calculating the auxiliary function.
module;
#include <cassert>
#include <stdexcept>
#include <string>

module qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;

static double Factorial[]= {1,1,2,6,24,120,720};

namespace qchem::BasisSet::Gaussian::Evaluators::Internal::Boys
{

AuxillaryFJ::AuxillaryFJ(int jMax) : itsjMax(jMax)
{
    // The Taylor step reads columns jmax..jmax+6 of the [121][23] table, so jmax+6 <= 22 is the real bound.
    static_assert(sizeof(theLookUp[0])/sizeof(double)==23, "Boys table width changed: re-derive thejMax");
    if (jMax<0 || jMax>thejMax)
        throw std::out_of_range("AuxillaryFJ: jmax="+std::to_string(jMax)+" exceeds the Boys F_j lookup table (jmax <= "
                               +std::to_string(thejMax)+", i.e. La+Lb+Lc+Ld <= "+std::to_string(thejMax)+")");
}

void AuxillaryFJ::GetFjAt(double T, rvec_t& Fj) const
{
    assert(Fj.size()==size_t(itsjMax+1)); // the size the constructor checked
    assert(T>=0                          );

    int jmax=itsjMax;
    double et=exp(-T);

    if(T<12)
    {
        int it=(int)(T*10.0+0.5);
        double Tstar=(double)it*0.1;
        double dt=Tstar-T;
        double* FT=&theLookUp[it][0];
        double fjmax=FT[jmax];
        for (int k=1; k<=6; k++) fjmax+=FT[k+jmax]*uintpow(dt,k)/Factorial[k];
        Fj[jmax]=fjmax;
        for (int j=jmax-1; j>=0; j--) Fj[j]=(2*T*Fj[j+1]+et)/(2*j+1);
    }
    if(T>=12 && T < 30.0)
    {
        double g=0;
        if(         T<15) g = 0.4999489092 - 0.2473631686/T + 0.321180909/(T*T) - 0.3811559346/(T*T*T);
        if(T>=15 && T<18) g = 0.4998436875 - 0.24249438  /T + 0.24642845 /(T*T);
        if(T>=18 && T<24) g = 0.499093162  - 0.2152832   /T;
        if(T>=24        ) g = 0.49;
        assert(g!=0.0); // g never got assigned.
        Fj[0] = 0.5*sqrt(Pi/T) - et*g/T;
        for (int j=0; j<jmax; j++) Fj[j+1]=((2*j+1)*Fj[j]-et)/(2*T);
    }
    if (T>=30.0)
    {
        Fj[0]=0.5*sqrt(Pi/T);
        for (int j=0; j<jmax; j++) Fj[j+1]=((2*j+1)*Fj[j]-et)/(2*T);
    }
}

} //namespace

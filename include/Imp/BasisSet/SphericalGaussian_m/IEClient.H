// File: SphericalGaussian_m/IEClient.H
#ifndef _SphericalGaussian_m_IEClient_H_
#define _SphericalGaussian_m_IEClient_H_

#include "Imp/BasisSet/SphericalGaussian/IEClient.H" 
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/BasisSet/BFGrouper.H"
#include "oml/vector.h"

namespace SphericalGaussian_m
{
    
typedef SphericalGaussian::IrrepIEClient IrrepIEClient;

struct IEData
{
    IEData() {}; 
    IEData(size_t N) : Ns(N), Ls(N), Ms(N), es(N), ns(N) {};
    Vector<size_t> Ns;  //Angular momentum.
    Vector<size_t> Ls;  //Angular momentum.
    Vector<int   > Ms;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};


//Use the same struct for a full basis set.;
struct IEClient 
    : public virtual ::IEClient 
    , public AtomIEClient
{
    void Append(const IrrepIEClient*);
    
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    
    virtual const SphericalGaussianCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SphericalGaussianCD*>(c);
    }

    virtual const IrrepIEClient* operator[](size_t i) const {return itsIrreps[i-1];}
    std::vector<const IrrepIEClient*> itsIrreps;

    
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_

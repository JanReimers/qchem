// File: PolarizedGaussian/IEClient.H
#ifndef _PolarizedGaussianIEClient_H_
#define _PolarizedGaussianIEClient_H_

#include "Polarization.H"
#include <IEClient.H> // <> prevents file search on this SphericalGaussian folder.
#include "oml/vector.h"
#include <vector>
#include <map>

namespace PolarizedGaussian
{

class Block;
class RadialFunction;


struct IEData
{
    void Init(std::vector<const Block*>&);
    std::vector<const RadialFunction*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    Vector<double>                     ns;      //Norm constants

};

struct IrrepIEClient
    : public virtual ::IrrepIEClient
    , public IEData
{
    IrrepIEClient() {};
    
    virtual size_t size() const {return radials.size();}
};

struct IEClient
    : public virtual ::IEClient
    , public IEData
{
    IEClient() {};
    void Append(const IrrepIEClient*);
    
    virtual size_t size() const {return radials.size();}
    bool SameIrrep(index_t ia, index_t ib) const;
    
    std::map <index_t,index_t> irrepIDs;  //<index,irrep#>
};
    
} //namespace PolarizedGaussian


#endif // _PolarizedGaussianIEClient_H_

// File: PolarizedGaussian/IEClient.H
#ifndef _PolarizedGaussianIEClient_H_
#define _PolarizedGaussianIEClient_H_

#include "Polarization.H"
#include <IEClient.H> // <> prevents file search on this SphericalGaussian folder.
#include "oml/vector.h"
#include <vector>

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
    void Init(std::vector<const Block*>& bs)
    {
        IEData::Init(bs);
        for (auto b:bs) blocks.push_back(b);
    }
    
    virtual size_t size() const {return radials.size();}
    std::vector<const Block*> blocks;
};
    
} //namespace PolarizedGaussian


#endif // _PolarizedGaussianIEClient_H_

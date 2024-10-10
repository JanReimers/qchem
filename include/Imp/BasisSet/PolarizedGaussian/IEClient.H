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

struct IrrepIEClient
    : public virtual ::IrrepIEClient
{
    
    IrrepIEClient() {};
    void Init(std::vector<const Block*>&);
    
    virtual size_t size() const {return radials.size();}

    
    std::vector<const RadialFunction*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    Vector<double>                     ns;      //Norm constants

};
    
} //namespace PolarizedGaussian


#endif // _PolarizedGaussianIEClient_H_

// File: Atom/l/Gaussian_IEC.H
#ifndef _Atom_l_Gaussian_IEC_H_
#define _Atom_l_Gaussian_IEC_H_

#include "../../AtomIEClient.H"

namespace Atoml
{
namespace Gaussian
{
    
    
struct IrrepIEClient
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N) : AtomIrrepIEClient(N) {};
private:
    virtual double Norm(double e, size_t l) const;
};

} //namespace
} //namespace

#endif // _Atom_l_Gaussian_IEC_H_

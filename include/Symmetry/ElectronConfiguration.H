// File: ElectronConfiguration.H
#ifndef _ElectronConfiguration_H_
#define _ElectronConfiguration_H_

#include <Symmetry/fwd.H>

class ElectronConfiguration
{
public:
    virtual int    GetN(const Irrep_QNs&) const=0;
    virtual void   Display() const=0;
};

#endif // _ElectronConfiguration_H_

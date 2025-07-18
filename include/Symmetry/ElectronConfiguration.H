// File: ElectronConfiguration.H
#ifndef _ElectronConfiguration_H_
#define _ElectronConfiguration_H_

import qchem.Symmetry.Irrep;

class ElectronConfiguration
{
public:
    virtual int    GetN(const Irrep_QNs&) const=0;
    virtual void   Display() const=0;
};

#endif // _ElectronConfiguration_H_

// File: Atom/radial/Slater/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom Slater Basis Sets.
#ifndef _Atom_radial_Slater_BS_Common_H_
#define _Atom_radial_Slater_BS_Common_H_

import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Common;

namespace Slater
{

// Common base handles all the radial aspects.
class BS_Common
: public ::BS_Common
, public ::AtomIE_BS_2E<double>
{
protected:
    virtual void Insert(bs_t* bs);
private:
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
};

}

#endif //_Atom_radial_Slater_BS_Common_H_
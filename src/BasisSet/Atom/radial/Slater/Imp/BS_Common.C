// File: Atom/radial/Slater/BS_Common.C  l/ml/kappa/mj independent part of BasisSet for atom Slater basis functions.
module;
#include <vector>
#include <iostream>
#include <cassert>
module qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Rk;
import qchem.IrrepBasisSet;
import qchem.BasisSet.Internal.IEClient;
import oml;

namespace Slater
{
    void BS_Common::Insert(bs_t* bs)
    {
        ::BS_Common::Insert(bs);
        auto iec=dynamic_cast<const IrrepIEClient*>(bs);
        assert(iec);
        IBS_Evaluator* eval=dynamic_cast<IBS_Evaluator*>(bs);
        assert(eval);
        Append(iec,eval);
    }
    
    
    
}
// File: BasisSet/Atom/Imp/BS.C Common for all atom basis sets.
module;
#include <cassert>
module qchem.BasisSet.Atom.BS;
import qchem.BasisSet.Internal.IEClient;

namespace Atom
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

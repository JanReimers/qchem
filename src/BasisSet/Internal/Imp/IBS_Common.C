// File: Imp/IBS_Common.C  Irrep Basis set common implementation.
module;
#include <iosfwd>
#include <iomanip>
#include <memory>
#include <cassert>
module qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisFunction;
import qchem.Symmetry;
import qchem.stl_io;

//-----------------------------------------------------------------------------
//
//  Construction zone
//

IBS_Common1::IBS_Common1(Symmetry* theQN)
    : itsSymmetry(theQN)
{
    assert(itsSymmetry);
};


//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
void IBS_Common1::Insert(bf_t* bf)
{
    assert(bf);
    itsBasisFunctions.push_back(std::shared_ptr<bf_t>(bf));
}







// File: BasisSet/Atom/BS.C Common for all atom basis sets.
export module qchem.BasisSet.Atom.BS;
import qchem.BasisSet.Internal.Common;
import qchem.BasisSet.Atom.IE;
import BasisSet.Atom.BS_Evaluator;

export namespace Atom
{
// Creates the Rk tool for HF ERIs
class BS_Common
: public ::BS_Common
, public ::AtomIE_BS_2E<double>
{
protected:
    BS_Common(BS_Evaluator* bse) : AtomIE_BS_2E<double>(bse) {};
    virtual void Insert(bs_t* bs);
};

}
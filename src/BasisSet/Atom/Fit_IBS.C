// File: BasisSet/Atom/Fit_IBS.C  Interface for a fitting Basis Set.
export module qchem.BasisSet.Atom.Fit_IBS;
export import qchem.Fit_IBS;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.IrrepBasisSet;

export namespace Atom
{
class Fit_IBS
    : public virtual ::Fit_IBS 
    , public Fit_IBS_Common
    , public AtomIE_Overlap<double>
    , public AtomIE_Fit
{
protected:
    Fit_IBS(const DB_cache<double>* db,const IE_Primatives* pie) 
    : AtomIE_Overlap<double>(db,pie)
    , AtomIE_Fit(db,pie)
    {};
};

} //namespace
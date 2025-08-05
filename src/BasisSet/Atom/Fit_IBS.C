// File: BasisSet/Atom/Fit_IBS.C  Interface for a fitting Basis Set.
export module qchem.BasisSet.Atom.Fit_IBS;
import qchem.Fit_IBS;
import qchem.BasisSet.Atom.IE;

export namespace Atom
{
class Fit_IBS
    : public virtual ::Fit_IBS 
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
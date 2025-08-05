// File: BasisSet/Atom/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron atom calculation.
export module qchem.BasisSet.Atom.Orbital_1E_IBS;
import qchem.Orbital_1E_IBS;
import qchem.BasisSet.Atom.IE;

export namespace Atom
{
template <class T> class Orbital_IBS
    : public virtual ::Orbital_IBS<T> //brings in Integrals_Overlap<T>
    , public AtomIE_Overlap<double>
    , public AtomIE_Kinetic<double>
    , public AtomIE_Nuclear<double>
{
protected:
    Orbital_IBS(const DB_cache<double>* db,const IE_Primatives* pie) 
    : AtomIE_Overlap<double>(db,pie)
    , AtomIE_Kinetic<double>(db,pie)
    , AtomIE_Nuclear<double>(db,pie) 
    {};
};

} //namespace
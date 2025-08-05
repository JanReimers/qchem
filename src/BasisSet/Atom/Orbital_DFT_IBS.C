// File: BasisSet/Orbital_DFT_IBS.C  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.


export module qchem.BasisSet.Atom.Orbital_DFT_IBS;
import qchem.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Atom.IE;

export namespace Atom
{
template <class T> class Orbital_DFT_IBS
    : public virtual ::Orbital_DFT_IBS<T> 
    , public Orbital_DFT_IBS_Common<double>
    , public AtomIE_DFT<double>
{
protected:
    Orbital_DFT_IBS(const DB_cache<double>* db,const IE_Primatives* pie) 
    : AtomIE_DFT<double>(db,pie)
    {};
};

} //namespace
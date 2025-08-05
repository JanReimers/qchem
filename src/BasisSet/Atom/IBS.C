// File: BasisSet/Atom/IBS.C Atom specific irrep basis sets.
export module qchem.BasisSet.Atom.IBS;
export import qchem.Orbital_1E_IBS;
export import qchem.Orbital_DFT_IBS;
export import qchem.Orbital_HF_IBS;
export import qchem.Fit_IBS;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.IrrepBasisSet;

export namespace Atom
{
template <class T> class Orbital_IBS
    : public virtual ::Orbital_IBS<T> //brings in Integrals_Overlap<T>
    , public Orbital_IBS_Common<double>
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

template <class T> class Orbital_HF_IBS
    : public virtual ::Orbital_HF_IBS<T> 
    , public Orbital_HF_IBS_Common<T>
{
protected:
    Orbital_HF_IBS(const DB_BS_2E<double>* db) : Orbital_HF_IBS_Common<T>(db) {};
};

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
// File: BasisSet/Atom/Orbital_HF_IBS.C  Interface for a Atom Hartree-Fock (HF) Orbital Irrep Basis Set.
export module qchem.BasisSet.Atom.Orbital_HF_IBS;
import qchem.Orbital_HF_IBS;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IrrepBasisSet;

export namespace Atom
{
template <class T> class Orbital_HF_IBS
    : public virtual ::Orbital_HF_IBS<T> 
    , public Orbital_HF_IBS_Common<T>
{
protected:
    Orbital_HF_IBS(const DB_BS_2E<double>* db) : Orbital_HF_IBS_Common<T>(db) {};
};

} //namespace
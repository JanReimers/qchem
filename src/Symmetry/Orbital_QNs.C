// File: Orbital_QNs.C  Encapsulate and sort a group of QNs associated with Atom/Molecule/LAttice orbitals.
module;
#include <iosfwd>

export module qchem.Symmetry.Orbital;
export import qchem.Symmetry.Irrep;
export import qchem.Types;

export struct Orbital_QNs
    : public Irrep_QNs
{   
    Orbital_QNs() : Irrep_QNs(), n(0) {};
    Orbital_QNs(size_t n, Spin ms,const sym_t& sym);
    Orbital_QNs(size_t n, const Irrep_QNs&);
    ~Orbital_QNs();

    virtual size_t  SequenceIndex() const; //Used for op<
    //size_t GetDegeneracy() const; //Same as Irrep

    std::ostream& Write(std::ostream&) const;

    size_t  n; //principle QN for this symmetry.
private:
    static const size_t n_max; //Used for calculating sequenxe indexes.
  
};
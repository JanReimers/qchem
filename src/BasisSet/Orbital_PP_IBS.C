// File: BasisSet/Orbital_PP_IBS.C  The species-field INTEGRAL SERVICE of an orbital basis.
module;
#include <map>
export module qchem.BasisSet.Orbital_PP_IBS;
export import qchem.BasisSet.SpeciesField;   // the argument vocabulary: SpeciesRadialField / SpeciesProjectorSet
import qchem.Structure;                       // Structure -- the positions and species the lattice sums run over
import qchem.Types;                           // hmat_t<T>

export namespace qchem::BasisSet
{

//! \brief The two integrals a species-attached external field asks of an orbital basis (V1.2; formerly
//! \c Pseudopotential::Integrals_Pseudo, whose arguments were pseudopotential types).  \tparam T the matrix
//! scalar.  A periodic basis assembles \f$\frac1\Omega\sum_a f(Z_a)\,e^{-i\Delta G\cdot\tau_a}\f$ from the
//! field's \f$q\f$ view; a real-space basis quadratures the \f$r\f$ view -- the SAME face, the basis owning HOW.
//! Nothing here is pseudopotential-specific: the arguments are a radial field and a projector set.
template <class T> class Orbital_PP_IBS
{
public:
    virtual ~Orbital_PP_IBS() {}
    //! \f$\langle i|\sum_a v_{Z_a}(|r-R_a|)|j\rangle\f$ for the given RANGE of \a field (Hermitian).  The
    //! \f$\Delta G=0\f$ term of a periodic assembly is the dropped uniform shift (\c CellMeanQ aligns it).
    virtual hmat_t<T> MakeSpeciesFieldMatrix(const Structure*, const SpeciesRadialField& field, FieldRange) const=0;
    //! \f$\langle i|\sum_a\sum_p w_p|\beta^a_p\rangle\langle\beta^a_p|\,|j\rangle\f$ (Hermitian), the \f$m\f$
    //! sum done by the basis (a plane-wave basis via the \f$(2l+1)P_l(\cos\gamma)\f$ addition theorem).
    virtual hmat_t<T> MakeProjectorMatrix(const Structure*, const SpeciesProjectorSet&) const=0;
    //! \brief The PER-ANGULAR-CHANNEL decomposition of \c MakeProjectorMatrix: one matrix per \f$l\f$, summing
    //! exactly to the whole.  DIAGNOSTIC (doc/SphericalLatticePlan.md I0: the MnO ordering campaign reads
    //! \f$E^{(l)}=\mathrm{Tr}(DV^{(l)})\f$ per channel).  Default: the whole lumped under \f$l=-1\f$, so a basis
    //! that has not implemented the split still answers correctly (callers treat \f$-1\f$ as "unresolved").
    virtual std::map<int,hmat_t<T>> MakeProjectorMatrixByL(const Structure* st, const SpeciesProjectorSet& ps) const
    {
        std::map<int,hmat_t<T>> byL;
        byL.emplace(-1, MakeProjectorMatrix(st, ps));
        return byL;
    }
};

} //namespace

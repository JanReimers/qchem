// File: Hamiltonian/Internal/PWTerms.C  Plane-wave (dcmplx) Kohn-Sham Hamiltonian terms.
//
// These are the THIN terms that complete the dependency inversion: each derives from the dcmplx term
// base (cStatic_HT/cDynamic_HT in qcHamiltonian), holds the abstract orbital basis cobs_t, dynamic_casts
// it UP to the abstract BasisSet::DFTPotential_IBS<dcmplx> capability (in qcBasisSet), and asks that high-
// level question -- "the external matrix", "the Hartree matrix for this density".  The basis owns the
// integration; the term owns no G-vectors or mesh.  Energies delegate to the density's DM_Contract.
module;
#include <iosfwd>
#include <memory>
export module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Hamiltonian.Internal.Term;        // cStatic_HT / cDynamic_HT + their _Imp cache bases
import qchem.BasisSet.DFTPotential_IBS;         // the abstract capability the terms dynamic_cast to
import qchem.Hamiltonian.Types;                 // cobs_t
import qchem.Structure;

export namespace qchem::Hamiltonian
{

//! External (pseudo)potential for a plane-wave basis (static, density-independent).  The pseudopotential
//! itself lives on the basis (set via PlaneWave_IBS::SetPseudopotential); this term just asks for the
//! assembled matrix.  Pair with the kinetic, Hartree and XC terms for a full Kohn-Sham Hamiltonian.
class PW_External
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> cl_t;
    PW_External(const cl_t& cl);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalculateMatrix(const cobs_t*, const Spin&) const;
    cl_t theStructure;
};

} //namespace

// File: BasisSet/Lattice_3D/APW_IBS.C  Augmented Plane Wave (APW) basis for one k-point.
//
// Lineage B, first IBS (see doc/PlaneWavePlan.md): plane waves in the interstitial region, matched
// (value-continuous) to free-particle radial solutions inside a single muffin-tin sphere at the cell
// origin.  This is the fixed-energy / empty-lattice DEMONSTRATOR: the APW secular matrix
//     Gamma(E) = H(E) - E O(E)
// is built per energy parameter E and is SINGULAR exactly at the free-electron energies E = 1/2|k+G|^2
// (the APW self-consistency condition).  Both H and O are energy-dependent (APW's nonlinearity), so
// rather than a single diagonalisation we expose Gamma(E) and test its singularity at the known levels.
//
// Derivation (empty lattice, V=0, symmetric-gradient variational form 1/2 int grad.phi* grad.phi'):
//   Gamma_{GG'} = (1/2 K_G.K_G' - E) [ delta_{GG'} - (4 pi R^2/Omega) j_1(|dG|R)/|dG| ]            (interstitial)
//               + (2 pi R^2/Omega) Sum_l (2l+1) P_l(cos g) j_l(K_G R) j_l(K_G' R) u_l'(R)/u_l(R)     (sphere)
// where K = k+G, cos g = K_G.K_G'/(|K_G||K_G'|), and for the empty lattice u_l(r)=j_l(qr), q=sqrt(2E),
// so u_l'(R)/u_l(R) = q j_l'(qR)/j_l(qR).  The radial integral collapses to this surface log-derivative
// (1/2 T_l - E R_l = 1/2 R^2 u_l(R) u_l'(R)), so no radial quadrature is needed.
//
// Framework note: APW_IBS shares the IrrepBasisSet<dcmplx> base with PlaneWave_IBS, but does NOT
// implement Orbital_1E_IBS -- its overlap is energy-dependent, so the (E-independent) 1E-integral
// interface does not apply.  Lineage B is a sequence of such IBS types (APW, LAPW, FLAPW, ...), each
// with its own interface needs.
module;
#include <iosfwd>
#include <string>
#include <vector>

export module qchem.BasisSet.Lattice_3D.APW_IBS;
import qchem.BasisSet.IrrepBasisSet;                 // IrrepBasisSet<dcmplx>
import qchem.BasisSet.Internal.IrrepBasisSetImp;     // GetSymmetry/GetSymt/GetIrrep
export import qchem.ReciprocalLattice;
import qchem.Types;

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief Augmented Plane Wave basis for one k-point (single muffin-tin sphere at the cell origin).
class APW_IBS
    : public virtual BasisSet::IrrepBasisSet<dcmplx>
    , public         BasisSet::IrrepBasisSetImp<dcmplx>
{
public:
    //! \param Ecut plane-wave cutoff for the interstitial G-set; \param Rmt muffin-tin radius (Bohr);
    //! \param lmax angular-momentum cutoff for the augmentation.
    APW_IBS(const ReciprocalLattice& recip, const ivec3_t& N, const ivec3_t& kIndex,
            double Ecut, double Rmt, size_t lmax);

    virtual size_t GetNumFunctions() const {return itsG.size();}

    //! \brief The APW secular matrix \f$\Gamma(E)=H(E)-E\,O(E)\f$ at energy parameter \a E (>0).
    //! Hermitian (real for the single origin sphere).  \f$\Gamma(E)\f$ is singular iff \a E is an APW
    //! eigenvalue; for the empty lattice those are the free-electron energies \f$\tfrac12|k+G|^2\f$.
    chmat_t MakeSecular(double E) const;

    // VectorFunction: interstitial plane-wave value.  The in-sphere augmentation is not evaluated
    // pointwise here (it needs Y_lm); the secular matrix does not require it.
    virtual cvec_t     operator() (const rvec3_t& r) const;
    virtual cvec3vec_t Gradient   (const rvec3_t& r) const;

    virtual std::string Name      () const {return "APW";}
    virtual std::string BasisSetID() const;
    virtual std::ostream& Write(std::ostream&) const;

private:
    ReciprocalLattice    itsRecip;  //!< Reciprocal cell (matrix B).
    rvec3_t              itsk;      //!< Fractional crystal momentum k = kIndex/N.
    double               itsVolume; //!< Direct cell volume Omega.
    double               itsRmt;    //!< Muffin-tin radius (Bohr).
    size_t               itsLmax;   //!< Angular-momentum cutoff.
    std::vector<ivec3_t> itsG;      //!< Reciprocal index triples within the cutoff.
};

} //namespace

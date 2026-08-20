// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/PGData.C  Flattened rep of PG IBS suitable for integral evaluation.
module;
#include <vector>
#include <string>

export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.Types;

export namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

struct PGData
{
    std::string BasisSetID() const; // geometry-aware cache identity: radial @ centre : pol per fn

    void Init(std::vector<const Block*>&);
    std::vector<const GaussianRF*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    rvec_t                             ns;      //Norm constants

    size_t size() const {return radials.size();}

    //! \brief Per-function SCREENING RADIUS: beyond it \f$|\chi_i(r)|<\varepsilon\f$ everywhere.
    //!
    //! A magnitude screen -- the only truncation this project allows -- expressed as the radius at which
    //! that magnitude bound is met, exactly as \c BetaSupportRadius does for the KB projectors.  The bound
    //! is \f$n_i\,d^{L}\,|R(d)|\f$: the contracted radial is spherically symmetric, and the Cartesian
    //! monomial obeys \f$|x^ay^bz^c|\le d^{L}\f$ over every direction, so this is an upper bound on the
    //! component at distance \a d and the radius is direction-independent.  Found by SCANNING (the radial
    //! is a sum of Gaussians and need not be monotone under contraction), taking the LAST \a d above
    //! \f$\varepsilon\f$ -- never the first below it.
    //!
    //! WHY IT EXISTS: this is the pointwise sweep behind every mesh quadrature (the atom-centred XC mesh's
    //! Phi tables, molecular Becke grids), and under a PERIODIC Bloch sum it runs once per LATTICE IMAGE
    //! per mesh point -- hundreds of image evaluations of which only the few nearby ones carry anything.
    //! Unscreened, the sweep pays a contracted exp() per shell for every one of them.
    //! Lazily built and cached (geometry-fixed, like the radials themselves).
    //! Built EAGERLY by \c Init (never lazily): the pointwise sweep runs inside \c MakePhi 's OpenMP
    //! region when \c GPW_OMP_THREADS is set, so a lazy first call would arrive on several threads at
    //! once and race on the resize.  Construction is single-threaded, which settles it with no flag --
    //! and a \c std::once_flag would also make \c PGData non-assignable (PG_LibCint assigns one).
    const rvec_t& Reaches() const {return itsReach;}

private:
    void   BuildReaches();
    rvec_t itsReach;

};

}
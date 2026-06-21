// File: BasisSet/Molecule/Evaluators/Internal/MnD/RNLM.C  Manage the R_NLM Hermite Coulomb auxiliaries.
//
// Generic MnD: addressed by a plain Hermite index (Index3), NOT the Cartesian Polarization -- that is the
// seam that lets PG_Cart_MnD, a future PG_Spherical_MnD, etc. share this.  (Polarization is implicitly an
// Index3, so Cartesian callers pass one unchanged.)
module;
#include <iosfwd>

export module qchem.BasisSet.Molecule.Evaluators.Internal.MnD.RNLM;
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD.Triangle3D;
import qchem.BasisSet.Molecule.Evaluators.Internal.MnD.Index3;
import qchem.Types;

export namespace BasisSet::Molecule::Evaluators::Internal::MnD
{

class RNLM
{
public:
    RNLM();
    // M&D eq 3.11 α=αₚ, dR = CP^2 for nuclear
    // M&D eq 3.32 α=αₚ*αQ/(αₚ+αQ) to electron repulsion
    //             dR = PQ^2
    RNLM(int Max, double α, const rvec3_t& dR);
    bool CheckLMax(int L) const {return L<=itsLMax;}

    double operator()(const Index3& p) const
    {
        return itsData(p.n,p.l,p.m);
    }

    // Element by element addition to this RNLM.
    void   Add  (const RNLM&, double Scale);
    void   Clear(                         )
    {
        itsData.Clear();
    }


private:
    int          itsLMax;
    Triangle3D   itsData;
};

} //namespace


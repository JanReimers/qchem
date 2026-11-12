// File: BasisSet/Gaussian/Evaluators/Internal/MnD/RNLM.C  Manage the R_NLM Hermite Coulomb auxiliaries.
//
// Generic MnD: addressed by a plain Hermite index (Index3), NOT the Cartesian Polarization -- that is the
// seam that lets PG_Cart_MnD, a future PG_Spherical_MnD, etc. share this.  (Polarization is implicitly an
// Index3, so Cartesian callers pass one unchanged.)
module;
#include <iosfwd>

export module qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM;
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D;
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Index3;
import qchem.Types;

export namespace qchem::BasisSet::Gaussian::Evaluators::Internal::MnD
{

class RNLM
{
public:
    RNLM();
    // M&D eq 3.11 α=αₚ, dR = CP^2 for nuclear
    // M&D eq 3.32 α=αₚ*αQ/(αₚ+αQ) to electron repulsion
    //             dR = PQ^2
    RNLM(int Max, double α, const rvec3_t& dR);
    //! Re-fill this RNLM in place for (Max, α, dR): the constructor's arithmetic exactly, reusing this object's storage.  A kernel that
    //! visits each primitive quartet once (the cold ERI4 build: 99.4 % of the global RNLM cache's lookups missed, OpenWork §4) holds one
    //! scratch RNLM and rebuilds it per quartet, with no cache descent, no insert and no heap traffic.
    void Build(int Max, double α, const rvec3_t& dR);
    bool CheckLMax(int L) const {return L<=itsLMax;}

    double operator()(const Index3& p) const
    {
        return itsData(p.n,p.l,p.m);
    }

    //! The Hermite-Coulomb triangle itself, for a kernel that walks its layout (see Triangle3D::Data).
    const Triangle3D& Table() const {return itsData;}

    // Element by element addition to this RNLM.
    void   Add  (const RNLM&, double Scale);
    void   Clear(                         )
    {
        itsData.Clear();
    }

    //! Heap bytes owned by the Hermite-Coulomb auxiliary store -- for the Cache2 RAM accounting.
    size_t HeapBytes() const {return itsData.HeapBytes();}


private:
    int          itsLMax;
    Triangle3D   itsData;
};

} //namespace


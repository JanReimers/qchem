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
import qchem.IntPow;   // uintpow / intpow (the recursion's base layer)

export namespace qchem::BasisSet::Gaussian::Evaluators::Internal::MnD
{
//! \brief THE HERMITE-COULOMB RECURSION (McMurchie-Davidson 4.1-4.4), ONE definition: from F_j(T), j = 0..Max, fill the layers
//! \f$R^{(j)}_{NLM}\f$, N+L+M <= Max-j, top layer first, through \a R(j,N,L,M) -> double&.  \c RNLM::Build runs it on its Triangle3D
//! scratch; the cache-free nuclear kernel Enzyme differentiates (GaussianRF::NuclearGradient) runs it on a flat local store.  INLINE
//! (a module's out-of-line body is not importable into the AD unit).  The base layer's powers of -2α are constants of the geometry.
template <class Fj_t, class R_t>
inline void HermiteCoulomb(int Max, double α, const rvec3_t& dR, const Fj_t& Fj, R_t&& R)
{
    R(Max,0,0,0)=uintpow(-2.0*α,Max)*Fj[Max];
    for (int j=Max-1; j>=0; j--)
    {
        R(j,0,0,0)=intpow(-2.0*α,j)*Fj[j];
        int n=Max-j;
        for (int M=0; M<n; M++)
        {
            R(j,0,0,M+1)= M>0 ? dR.z*R(j+1,0,0,M)+M*R(j+1,0,0,M-1) : dR.z*R(j+1,0,0,M);
            for (int L=0; L<n-M; L++)
            {
                R(j,0,L+1,M)= L>0 ? dR.y*R(j+1,0,L,M)+L*R(j+1,0,L-1,M) : dR.y*R(j+1,0,L,M);
                for (int N=0; N<n-M-L; N++)
                {
                    R(j,N+1,L,M)= N>0 ? dR.x*R(j+1,N,L,M)+N*R(j+1,N-1,L,M) : dR.x*R(j+1,N,L,M);
                }
            }
        }
    }
}


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


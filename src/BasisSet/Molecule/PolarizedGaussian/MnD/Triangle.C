// File: BasisSet/Molecule/PolarizedGaussian/MnD/Triangle3D.C  A 3-index triangular data structure.
module;
#include <iosfwd>
#include <vector>
export module qchem.BasisSet.Molecule.PolarizedGaussian.MnD.Triangle3D;
import qchem.Streamable;

#if DEBUG
#define CHECK(i,j,k) Check(i,j,k)
#else
#define CHECK(i,j,k)
#endif
//
//  Implements an triangle data structure, where N+L+M is always < MaxSum.
//
export class Triangle3D
    : public virtual Streamable
{
public:
    Triangle3D(   );
    Triangle3D(int);
    Triangle3D& operator=(const Triangle3D&);

    double operator()(int i,int j, int k) const
    {
        CHECK(i,j,k);
        return itsData[(N-i)*(N-i+1)*(N-i+2)/6+(N-i-j)*(N-i-j+1)/2+N-i-j-k];
    }
    double& operator()(int i,int j, int k)
    {
        CHECK(i,j,k);
        return itsData[(N-i)*(N-i+1)*(N-i+2)/6+(N-i-j)*(N-i-j+1)/2+N-i-j-k];
    }
    // Element by element addition of another triangle.
    void   Add  (const Triangle3D&, double Scale);
    void   Clear(                             )
    {
        itsData.clear();
        N=-1;
    }

    virtual std::ostream&  Write(std::ostream&) const;

private:
    void    Check(int,int,int) const;

    int     N;
    std::vector<double> itsData;
};

#undef CHECK


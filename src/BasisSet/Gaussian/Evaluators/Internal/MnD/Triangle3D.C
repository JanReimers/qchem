// File: BasisSet/Gaussian/Evaluators/Internal/MnD/Triangle3D.C  A 3-index triangular data structure.
module;
#include <iosfwd>
#include <vector>
export module qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D;
import qchem.Streamable;

#if DEBUG
#define CHECK(i,j,k) Check(i,j,k)
#else
#define CHECK(i,j,k)
#endif

export namespace qchem::BasisSet::Gaussian::Evaluators::Internal::MnD
{
//
//  Implements an triangle data structure, where N+L+M is always < MaxSum.
//
class Triangle3D
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
    //! RAW ACCESS for a hot loop that walks the layout itself (Repulsion4CTab hoists the offset arithmetic out of its Hermite loops,
    //! OpenWork §4): element (i,j,k) is Data()[Layer(N-i) + Row(N-i-j) + (N-i-j-k)] with N = MaxSum(), exactly operator()'s formula.
    const double* Data   () const {return itsData.data();}
    int           MaxSum () const {return N;}
    static constexpr size_t Layer(int x) {return size_t(x)*(x+1)*(x+2)/6;}
    static constexpr size_t Row  (int x) {return size_t(x)*(x+1)/2;}

    // Element by element addition of another triangle.
    void   Add  (const Triangle3D&, double Scale);
    void   Clear(                             )
    {
        itsData.clear();
        N=-1;
    }

    virtual std::ostream&  Write(std::ostream&) const;

    //! Heap bytes owned by the backing store -- for the RNLM/Cache2 RAM accounting.
    size_t HeapBytes() const {return itsData.size()*sizeof(double);}

private:
    void    Check(int,int,int) const;

    int     N;
    std::vector<double> itsData;
};

} //namespace

#undef CHECK


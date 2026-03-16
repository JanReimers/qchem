// File ERI4T.C  Generical containter for and ERI (Electron Repulsion Integral) 4 index super matrix.
module;
export module qchem.BasisSet.Internal.ERI4T;
import oml;
export import qchem.Types;
//------------------------------------------------------------------
//
// Two electron repulsion integral have the following form:
//          /
// <ij|kl>= | g_i(r1)*g_j(r1)*1/r_12 * g_k(r2)*g_l(r2) * d^3 r_2
//          /
//
// where g is a basis function. There is obviously some symmetry with regard to the
// i,j,k,l indicies.  In particular, the follwing swaps are allowed
//   i <--> j
//   k <--> l
//   i,j <--> k,l
// resulting in a eightfold storage reduction.
//

export template <class T,template<class> class M> class ERI4T
{
public:
    ERI4T() {};
    ERI4T(size_t Nab, size_t Ncd);
    const M<T>& operator()(size_t a, size_t b) const {return itsData.ref(a,b);}
          M<T>& operator()(size_t a, size_t b)       {return itsData(a,b);}
    

    size_t size() const;
    size_t Nab() const {return itsData.GetNumRows();}
    size_t Ncd() const {return itsData(1,1).GetNumRows();}
    MatLimits GetLimits() const {return itsData.GetLimits();}
    auto rows() const {return itsData.rows();}
    auto cols() const {return itsData.cols();}
    auto rows(size_t i) const {return itsData.rows(i);}
    auto cols(size_t i) const {return itsData.cols(i);}
private:
    M<M<T> > itsData;
};

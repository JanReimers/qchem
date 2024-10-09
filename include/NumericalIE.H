// File: IntegralEngine.H  Abtract integral engine interface.
#ifndef _IntegralEngine_H_
#define _IntegralEngine_H_


#include "types.H"


//--------------------------------------------------------------------------------
//
//  This class is responsible for calculating all integrals over basis functions.
//  The method of integral evaluation is of course strongly dependant on the
//  precise details of basis function or basis set.  Hence for each concrete type
//  of basis function there is a corresponding basis set and integral engine.
//  These classes are nessecarily tightly coupled.
//
//  Integrals supprted are:
//
//  Overlap                            <i|j>
//  Overlap  (BasisSet      &     b)   <i|bj>
//  Overlap  (ScalarFunction&     f)   <i|f>   *** f is NOT normalized ***
//  Overlap3C(BasisSet&           b)   <ij|b>  for all basis functions b, returns an optr_vector1<SMatrix>.
//
//  Repulsion                            <i(1)|1/r12|j(2)>
//  Repulsion  (BasisSet      &     b)   <i(1)|1/r12|bj(2)>
//  Repulsion  (ScalarFunction&     f)   <i(1)|1/r12|f(2)>      *** f is NOT normalized ***
//  Repulsion3C(BasisSet&           b)   <i(1)j(1)|1/r12|b(2)>  for all basis functions b, returns an optr_vector1<SMatrix>.
//  Repulsion4C(                     )   <i(1)j(1)|1/r12|k(2)l(2)>  for all internal basis functions.
//  Normalization        1/sqrt(<i|i>)
//  Charge               <i|1>
//
//  Kinetic              <i|-1/2grad^2|j>
//  Nuclear(Cluster& cl) sum over atoms c {<i|1/(r-Rc)|j>}
//

//
//  Set this define if you want use this engine for debugging analytic integrals.
//
//#define USE_FOR_DEBUGGING_ANALYTIC
//
//  only build functions the have UT coverage.
//
#define UT_COVERAGE_ONLY

template <class T> class NumericalIE
{
public:
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector <double> RVec;
    typedef SMatrix<double> RSMat;
    typedef const TIrrepBasisSet<T> bs_t;
    
    typedef const ScalarFunction<double> Rf;
    typedef const ScalarFunction<T>      Sf;
    typedef const VectorFunction<T>      Vf;

    virtual ~NumericalIE()  {};

    virtual RVec MakeNormalization(bs_t& a) const=0;
    virtual Vec  MakeOverlap    (Vf&    , Rf&) const=0;
    virtual Vec  MakeRepulsion  (Vf&    , Rf&) const=0;
    //
    //  These currently have no UT coverage.
    //
#ifndef UT_COVERAGE_ONLY
    virtual RVec MakeCharge       (bs_t& a) const=0;
    virtual SMat MakeOverlap      (bs_t& a) const=0;
    virtual Mat  MakeOverlap      (bs_t& a,bs_t& b) const=0;
    virtual Mat  MakeRepulsion    (bs_t& a,bs_t& b) const=0;
    virtual SMat MakeRepulsion    (bs_t& a) const=0;
#endif

    
    
#ifdef USE_FOR_DEBUGGING_ANALYTIC
    virtual void MakeOverlap3C  (ERI3&,bs_t&   ) const=0;
    virtual void MakeRepulsion3C(ERI3&,bs_t&   ) const=0;
    virtual void MakeRepulsion4C(ERI4& Coulomb, ERI4& exchange, const BasisGroup*) const=0;
    virtual SMat MakeKinetic(              ) const=0;
    virtual SMat MakeNuclear(const Cluster&) const=0;
#endif

    virtual NumericalIE* Clone  (         ) const=0;
//    static  NumericalIE* Factory(std::istream&);
};


#endif

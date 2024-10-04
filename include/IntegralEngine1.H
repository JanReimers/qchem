// File: IntegralEngine1.H  Abtract integral engine interface.
#ifndef _IntegralEngine1_H_
#define _IntegralEngine1_H_


#include "types.H"
#include "Misc/pmstream.h"
#include <vector>
#include "Misc/ERIList1.H"


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
template <class T> class IntegralEngine1
    : public virtual PMStreamableObject
{
public:
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector <double> RVec;
    typedef SMatrix<double> RSMat;
    typedef MatrixList<T>  MList;
    typedef std::vector<const IntegralEngine1<double>*> iev_t;
    typedef std::pair<ERIList1,ERIList1> jk_t;


    virtual ~IntegralEngine1()  {};
    
    virtual size_t size() const=0;

    virtual SMat MakeOverlap   (                             ) const=0;
    virtual Mat  MakeOverlap   (const IntegralEngine1*       ) const=0;
    virtual Vec  MakeOverlap   (const ScalarFunction<double>&) const=0;
    virtual void MakeOverlap3C  (MList&,const IntegralEngine1*   ) const=0;

    virtual SMat MakeRepulsion  (                             ) const=0;
    virtual Mat  MakeRepulsion  (const IntegralEngine1*       ) const=0;
    virtual Vec  MakeRepulsion  (const ScalarFunction<double>&) const=0;
    virtual void MakeRepulsion3C(MList&,const IntegralEngine1*   ) const=0;
    virtual void MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const iev_t&) const=0;
    virtual jk_t Make4C(const iev_t&) const=0;

    virtual RVec MakeNormalization() const=0;
    virtual RVec MakeCharge       () const=0;
//
    virtual SMat MakeKinetic(              ) const=0;
    virtual SMat MakeNuclear(const Cluster&) const=0;
//    virtual RSMat MakeInverse(const RSMat& ) const=0; //Numerically stable algo required.
//
    virtual std::ostream& Write(std::ostream&) const=0;
    virtual std::istream& Read (std::istream&)      =0;
    virtual IntegralEngine1* Clone  (         ) const=0;
    static  IntegralEngine1* Factory(std::istream&);
};


#endif

// File: AnalyticIE.H  Abtract integral engine interface.
#ifndef _AnalyticIE_H_
#define _AnalyticIE_H_


#include "types.H"
#include "Misc/pmstream.h"
#include <vector>
#include "Misc/ERI4.H"


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
template <class T> class AnalyticIE
    : public virtual PMStreamableObject
{
public:
    typedef Matrix<T>      Mat;
    typedef SMatrix<T>     SMat;
    typedef Vector<T>      Vec;
    typedef Vector <double> RVec;
    typedef SMatrix<double> RSMat;
    typedef std::vector<SMat> ERI3;
    typedef std::pair<ERI4,ERI4> jk_t;
    typedef const IrrepIEClient iec_t;
    typedef std::vector<iec_t*> iecv_t;


    virtual ~AnalyticIE()  {};
    
    virtual size_t size() const=0;

    virtual RVec MakeNormalization() const=0;
    virtual RVec MakeCharge       (iec_t* a) const=0;
    
    virtual  SMat MakeOverlap  (iec_t* a        ) const=0;
    virtual  SMat MakeKinetic  (iec_t* a        ) const=0;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const=0;
    virtual  SMat MakeRepulsion(iec_t* a         ) const=0; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(const AnalyticIE*,const AnalyticIE*) const=0; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t* a,iec_t* b) const=0; //Used for charge and Vxc fitting.
    static  RSMat MakeInverse  (const RSMat&     ); //Numerically stable algo required.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const=0; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const=0; //Used for DFT
    virtual jk_t Make4C(const iecv_t&) const=0; //Used for Hartree-Fock

    virtual std::ostream& Write(std::ostream&) const=0;
    virtual std::istream& Read (std::istream&)      =0;
    virtual AnalyticIE* Clone  (         ) const=0;
    static  AnalyticIE* Factory(std::istream&);
};


#endif //_AnalyticIE_H_

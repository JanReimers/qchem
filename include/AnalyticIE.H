// File: AnalyticIE.H  Abtract integral engine interface.
#ifndef _AnalyticIE_H_
#define _AnalyticIE_H_


#include "types.H"
#include <vector>
#include <map>
#include <iosfwd>

class IrrepIEClient;
//--------------------------------------------------------------------------------
//
//!  \brief This class is abstract interface used for calculating integrals over the basis functions.
//!
//!  The method of integral evaluation is of course strongly dependant on the
//!  precise details of basis function or basis set.  Hence for each concrete type
//!  of basis function there is a corresponding basis set and integral engine.
//!  The integral engine does not see the basis function/set objects directly. Instead it
//!  works with the IrrepIEClient interface which simply supply list of
//!  exponents, polarizations and normalization constants for the particular irrep basis set.
//!  All functions except MakeNormalization return normalized integrals.
//!
template <class T> class AnalyticIE
{
public:
    typedef Matrix <T>      Mat;
    typedef SMatrix<T>      SMat;
    typedef Vector <T>      Vec;
    typedef Vector <double> RVec;
    typedef SMatrix<double> RSMat;
    typedef std::vector<SMat> ERI3;
    typedef const IrrepIEClient iec_t;
    typedef std::map<size_t,std::map<size_t,ERI4> > erij_t;

    virtual ~AnalyticIE()  {};
    virtual void Append(const IrrepIEClient*);

    //! Normalization constants \f$ \frac{1}{\sqrt{\int d^{3}\vec{r}\:g_{a}^{2}\left(\vec{r}\right)}} \f$
    //virtual RVec MakeNormalization(const IrrepIEClient* a) const=0;
    //! Total charge \f$ \left\langle a\left|\right|1\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right) \f$
    virtual RVec MakeCharge       (const IrrepIEClient* a) const=0;
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    virtual  SMat MakeOverlap  (const IrrepIEClient* g) const=0;
    //! Kinetic energy \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual  SMat MakeKinetic  (const IrrepIEClient* g) const=0;
    //! Kinetic energy cross \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}f_{b}\left(\vec{r}\right)\f$
    virtual   Mat MakeKinetic  (const IrrepIEClient* g, const IrrepIEClient* f) const=0;
    //! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
    virtual  SMat MakeNuclear  (const IrrepIEClient* g,const Cluster&   ) const=0;
    //! Rest mass \f$ \left\langle a\left|\left(\beta-\alpha\right)c^{2}\right|b\right\rangle =\left(\beta-\alpha\right)c^{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$
    virtual  SMat MakeRestMass (const IrrepIEClient* g) const=0;
    //! Single basis set repulsion \f$ \left\langle a\left(\vec{r}_{1}\right)\left|\frac{1}{r_{12}}\right|b\left(\vec{r}_{2}\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{b}\left(\vec{r}_{2}\right)\f$
    virtual  SMat MakeRepulsion(const IrrepIEClient* g) const=0; //Used for charge and Vxc fitting.
    //! Cross basis set repulsion \f$ \left\langle a\left(\vec{r}_{1}\right)\left|\frac{1}{r_{12}}\right|b\left(\vec{r}_{2}\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{b}\left(\vec{r}_{2}\right) \f$
    virtual   Mat MakeRepulsion(const IrrepIEClient* g,const IrrepIEClient* f) const=0; //Used for charge and Vxc fitting.

    //! 3 centre overlap used for DFT \f$ \left\langle ab\left|1\right|c\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right)f_{c}\left(\vec{r}\right) \f$
    virtual ERI3 MakeOverlap3C  (const IrrepIEClient* g,const IrrepIEClient* f) const=0; //Used for DFT
    //! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
    virtual ERI3 MakeRepulsion3C(const IrrepIEClient* g,const IrrepIEClient* f) const=0; //Used for DFT
    //! \brief 4 center Coulomb and exchange used for Hartree-Fock.
    //! \f$ J_{ab,cd}=\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)d\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}g_{c}\left(\vec{r}_{2}\right)g_{d}\left(\vec{r}_{2}\right) \f$.
    //! \f$ K_{ab,cd}=\left\langle a\left(1\right)d\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)b\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{d}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}g_{c}\left(\vec{r}_{2}\right)g_{b}\left(\vec{r}_{2}\right) \f$. \n
    //! These integrals are not block diagonal for each irreducible representation. Therefore they
    //! need to be evaluated over a full basis set which is a list of IrrepIEClients.
    virtual void MakeDirect  (erij_t& Jab) const; //Used for Hartree-Fock
    virtual void MakeExchange(erij_t& Kab) const; //Used for Hartree-Fock
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const=0;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const=0;
    
    //! \brief Return the Penrose inverse of a symmetric matrix using SVD decomposition
    //! If \f$ S=UsV^{\dagger} \f$, then \f$ S^{-1}=V\frac{1}{s}U^{\dagger} \f$
    static  RSMat MakeInverse  (const RSMat&,const LAParams&); //Numerically stable algo required.
private:
    std::vector<const IrrepIEClient*> itsIrreps; //Used for 2-electron integrals.

};

//
// Interface for 1 electron integrals.  Can be used by Integral engines and Integral data bases.
// Two of these are only used for the Dirac Hamiltonian.  SHould that be a separate interface?
//
template <class T> class Integrals_1E
{
public:
    //! Instance then ensures integrals are calculated once and then buffered.
    static Integrals_1E* theDB; 

    typedef const TIrrepBasisSet<T> bs_t;
    typedef  Matrix<T>  Mat;
    typedef const SMatrix<T>& SMat_ref;
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    //! Kinetic energy \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    //! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
    //! Rest mass \f$ \left\langle a\left|\left(\beta-\alpha\right)c^{2}\right|b\right\rangle =\left(\beta-\alpha\right)c^{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$
    virtual SMat_ref Integrals(qchem::IType,const bs_t* g,const Cluster* cl=0) const=0;
       
};

#endif //_AnalyticIE_H_

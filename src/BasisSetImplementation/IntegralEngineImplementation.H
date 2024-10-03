// File: IntegralEngineImplementation.H  Implement common attributes for Integral Engines.
#ifndef _IntegralEngineImplementation_H_
#define _IntegralEngineImplementation_H_

#include "IntegralEngine.H"
#include "Misc/Persistent/IDRef.H"
#include "oml/vector.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <vector>

//--------------------------------------------------------------------------------
//
//  Normalization and total charge integrals are provided by the basis function
//  Interface, and can be implemented here in an abstract way.
//
template <class T> class IntegralEngineImplementation
    : virtual public IntegralEngine<T>
{
protected:
    typedef IntegralEngine<T> Base;
    typedef typename Base::Vec  Vec;
    typedef typename Base::Mat  Mat;
    typedef typename Base::RVec RVec;
    typedef typename Base::SMat SMat;
    typedef typename Base::RSMat RSMat;
public:
    IntegralEngineImplementation();
    IntegralEngineImplementation(const IntegralEngineImplementation&);
    virtual void Insert(const TBasisSet<T>*);

    virtual RVec MakeNormalization  () const;
    virtual RVec MakeCharge         () const;
    
    virtual RSMat MakeInverse(const RSMat&) const;
    
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void CheckBasisSet() const;         //Check that itsBasisSet has been initialized.
    void Normalize(SMat& m) const
    {
        auto nn=OuterProduct(itsNormalizations);
        m=DirectMultiply(m,nn);
//        ::Normalize(m,itsNormalizations);
    }
    void Normalize(Vec& v) const
    {
        v=DirectMultiply(v,itsNormalizations);
    }
    void Normalize(const RVec&, Mat&, const RVec&) const;
    std::vector<SMat> MakeMatrixList (int n) const;

    IDRef<const TBasisSet<T> > itsBasisSet;       //Permenant hook into the basis set.
    size_t                     itsN;              //Number of basis functions.
    RVec                       itsNormalizations; //Normalization constants.

private:
    void InitNormalizations();

};


#endif //_IntegralEngineImplementation_H_

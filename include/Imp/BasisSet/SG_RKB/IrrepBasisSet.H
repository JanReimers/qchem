// File: SphericalGaussian_RKB_IrrepBasisSet.H  .
#ifndef _SphericalGaussian_RKB_BS_H_
#define _SphericalGaussian_RKB_BS_H_



#include "IEClient.H"
#include "Imp/BasisSet/Dirac_IBS.H"
#include "Imp/BasisSet/SG_RKB/IntegralEngine.H"
#include "../IrrepCommon.H"
#include "../TIrrepCommon.H"
#include "../SphericalGaussian/IntegralEngine.H"
#include "oml/vector.h"
class DiracIntegralTests;


namespace SphericalGaussian_RKB
{

template <class T> class Large_IE
: public virtual Integrals_RKBL<T>
, public DB_RKBL<T>
{
public:
    typedef Integrals_RKBL<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::Mat  Mat;

    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const;
    virtual Mat  MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
};

template <class T> class Large_Orbital_IBS
    : public virtual ::Orbital_RKBL_IBS<T>
    , public     IrrepBasisSetCommon
    , protected  TIrrepBasisSetCommon<T> 
    , public     Large_IE<T>
    , public     IrrepIEClient
{
    public:
    Large_Orbital_IBS() {};
    Large_Orbital_IBS(const LAParams&,IntegralDataBase<T>*,
        const Vector<T>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};

template <class T> class Small_IE
: public virtual Integrals_RKBS<T>
, public DB_RKBS<T>
{
public:
    typedef Integrals_RKBL<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::Mat  Mat;
    virtual T    Integral(qchem::IType,double ea , double eb,size_t l) const;
};

template <class T> class Small_Orbital_IBS
    : public virtual ::Orbital_RKBS_IBS<T>
    , public     IrrepBasisSetCommon
    , protected  TIrrepBasisSetCommon<T> 
    , public     Small_IE<T>
    , public     Small_IrrepIEClient
{
public:
    Small_Orbital_IBS() {};
    Small_Orbital_IBS(const LAParams&,IntegralDataBase<T>*,const Large_Orbital_IBS<T>*);
    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};


class Large_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         Orbital_IBS_Common<double>
    , public         SphericalGaussian::IntegralEngine1
    , public         IrrepIEClient
{
public:
    Large_IrrepBasisSet();
    Large_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,
        const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}

    virtual IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;
};
class Small_IrrepBasisSet
    : public virtual TIrrepBasisSet<double>
    , public         IrrepBasisSetCommon
    , public         Orbital_IBS_Common<double>
    , public         Small_IrrepIEClient
    , public         SphericalGaussian_RKB::Small_IntegralEngine1

{
public:
    Small_IrrepBasisSet() {};
    Small_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Large_IrrepBasisSet*);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual ::IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

};

class Dirac_IrrepBasisSet
    : public virtual Orbital_IBS<double>
    , public Dirac::IrrepBasisSet<double>
    , public         Dirac_IrrepIEClient
{
    typedef DB_1E<double>::SMat SMat;
public:
    Dirac_IrrepBasisSet() {};
    Dirac_IrrepBasisSet(const LAParams&,IntegralDataBase<double>*,const Vector<double>& exponents, int kappa);

    virtual UniqueID::IDtype     GetID() const {return UniqueID::GetID();}
    virtual ::IrrepBasisSet* CreateCDFitBasisSet(const Cluster*) const;
    virtual ::IrrepBasisSet* CreateVxcFitBasisSet(const Cluster*) const;

    virtual std::ostream&  Write(std::ostream&    ) const;
    virtual std::istream&  Read (std::istream&    )      ;
    virtual ::IrrepBasisSet* Clone(            ) const;
    virtual ::IrrepBasisSet* Clone(const RVec3&) const;

private:
    friend class ::DiracIntegralTests;

};


} //namespace

#endif //_SphericalGaussian_RKB_BS_H_

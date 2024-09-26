// File: TOrbitalGroupImplementation.H  general orbital group implementation.
#ifndef _TOrbitalGroupImplementation_H_
#define _TOrbitalGroupImplementation_H_



#include "Orbital/TOrbitalGroup.H"
#include "OrbitalImplementation/OrbitalGroupImplementation.H"
#include "BasisSet/TBasisSet.H"

template <class T> class TOrbitalGroupBrowser;
template <class T> class TOrbitalGroupIterator;

template <class T> class TOrbitalGroupImplementation
    : public virtual TOrbitalGroup<T>
    , private OrbitalGroupImplementation
{
    typedef VectorFunction<T> Base;
    typedef typename Base::Mat     Mat;  //Matrix.
    typedef typename Base::SMat    SMat; //Symmetrix matrix.
    typedef typename Base::Vec     Vec;  //Vector of scalars.
    typedef typename Base::Vec3    Vec3;   //3 vector (possibly complex).
    typedef typename Base::Vec3Vec Vec3Vec;//vector of 3 space vectors.
    typedef typename Base::RVec3   RVec3;  //Real space vector.
    typedef typename Base::RVec    RVec;
    typedef typename Base::Vec3Mat Vec3Mat;//matrix of 3 space vectors.
public:
    TOrbitalGroupImplementation();
    TOrbitalGroupImplementation(const rc_ptr<const BasisSet>&,
                                const Mat&,
                                const RVec&,
                                const Spin&);
    virtual ~TOrbitalGroupImplementation() {};

    typedef Vector3D<T> TVec;

    virtual Vec     operator()(const RVec3&) const;
    virtual Vec3Vec Gradient  (const RVec3&) const;

    virtual ChargeDensity* GetChargeDensity(Spin) const;

    virtual std::ostream&          Write(std::ostream&) const;
    virtual std::istream&          Read (std::istream&)      ;

private:
    TOrbitalGroupImplementation(const TOrbitalGroupImplementation&);

    friend class TOrbitalGroupBrowser<T>;
    friend class TOrbitalGroupIterator<T>;

    virtual SMat CalculateDensityMatrix() const;
};


#endif //_TOrbitalGroupImplementation_H_

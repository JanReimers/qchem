// File: PlaneWaveIE.H  Integral Engine for polarized gaussians (contracted or primative).
#ifndef _PlaneWaveIE_H_
#define _PlaneWaveIE_H_



#include "BasisSetImplementation/IntegralEngineImplementation.H"
#include "Cluster/UnitCell.H"

class PlaneWaveBS;
template <class T> class SList;

class PlaneWaveIE
    : public virtual IntegralEngine<std::complex<double> >
    , private IntegralEngineImplementation<std::complex<double> >
{
public:
    PlaneWaveIE(                  );
    PlaneWaveIE(const PlaneWaveIE&);
    virtual void Insert(TBasisSet<std::complex<double> >*);

    virtual SMat MakeOverlap    (                                         ) const;
    virtual Mat  MakeOverlap    (const TBasisSet<std::complex<double> >&       ) const;
    virtual Vec  MakeOverlap    (const ScalarFunction<double>&            ) const;
    virtual void MakeOverlap3C  (MList&,const TBasisSet<std::complex<double> >&) const;

    virtual SMat MakeRepulsion  (                                         ) const;
    virtual Mat  MakeRepulsion  (const TBasisSet<std::complex<double> >&       ) const;
    virtual Vec  MakeRepulsion  (const ScalarFunction<double>&            ) const;
    virtual void MakeRepulsion3C(MList&,const TBasisSet<std::complex<double> >&) const;
    virtual void MakeRepulsion4C(ERIList& , ERIList&,const BasisGroup* ) const;

    virtual SMat MakeKinetic    (              ) const;
    virtual SMat MakeNuclear    (const Cluster&) const;

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual IntegralEngine<std::complex<double> >* Clone() const;

private:
    std::vector<RVec3> itsGs;
    RVec3              itsK;  //Bloch wave vector.
    UnitCell           itsRLCell;
};


#endif

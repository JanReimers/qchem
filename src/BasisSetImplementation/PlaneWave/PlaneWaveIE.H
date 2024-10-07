// File: PlaneWaveIE.H  Integral Engine for polarized gaussians (contracted or primative).
#ifndef _PlaneWaveIE_H_
#define _PlaneWaveIE_H_


#include "AnalyticIE.H"
#include "Cluster/UnitCell.H"

class PlaneWaveBS;
template <class T> class SList;

class PlaneWaveIE
    : public virtual AnalyticIE<std::complex<double> >
{
    typedef AnalyticIE<std::complex<double>> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;
    typedef IE::iev_t iev_t;
    typedef IE::jk_t jk_t;
public:
    PlaneWaveIE(                  );
    PlaneWaveIE(const PlaneWaveIE&);
    virtual void Insert(TBasisSet<std::complex<double> >*);

    virtual size_t size() const {return 0;}
    
    virtual RVec MakeNormalization() const;
    virtual RVec MakeCharge       () const;
    
    virtual  SMat MakeOverlap  (                 ) const;
    virtual  SMat MakeKinetic  (                 ) const;
    virtual  SMat MakeNuclear  (const Cluster&   ) const;
    virtual  SMat MakeRepulsion(                 ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(const AnalyticIE*) const; //Used for charge and Vxc fitting.

    virtual void MakeOverlap3C  (ERI3&,const AnalyticIE*) const; //Used for DFT
    virtual void MakeRepulsion3C(ERI3&,const AnalyticIE*) const; //Used for DFT
    virtual jk_t Make4C(const iev_t&) const; //Used for Hartree-Fock
    

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual AnalyticIE<std::complex<double> >* Clone() const;

private:
    std::vector<RVec3> itsGs;
    RVec3              itsK;  //Bloch wave vector.
    UnitCell           itsRLCell;
};


#endif

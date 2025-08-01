// File: PlaneWaveIE.H  Integral Engine for polarized gaussians (contracted or primative).
#ifndef _PlaneWaveIE_H_
#define _PlaneWaveIE_H_


import qchem.BasisSet.Internal.Integrals;

import Cluster.UnitCell;

class PlaneWaveBS;
template <class T> class SList;

class PlaneWaveIE
    : public virtual AnalyticIE<std::complex<double> >
{
    typedef AnalyticIE<std::complex<double>> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::iev_t iev_t;
    typedef IE::jk_t jk_t;
public:
    PlaneWaveIE(                  );
    PlaneWaveIE(const PlaneWaveIE&);
    virtual void Insert(TIrrepBasisSet<std::complex<double> >*);

    virtual size_t size() const {return 0;}
    
    virtual RVec MakeNormalization() const;
    virtual RVec MakeCharge       () const;
    
    virtual  SMatrix<T> MakeOverlap  (                 ) const;
    virtual  SMatrix<T> MakeGrad2  (                 ) const;
    virtual  SMatrix<T> MakeNuclear  (const Cluster&   ) const;
    virtual  SMatrix<T> MakeRepulsion(                 ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(const AnalyticIE*,const AnalyticIE*) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (const AnalyticIE*) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(const AnalyticIE*) const; //Used for DFT
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

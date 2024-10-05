// File: RNLM.H  Manager the RNLM Auxillary functions.
#ifndef _RNLM_H_
#define _RNLM_H_



#include "Misc/Polarization.H"
#include "Misc/Triangle/Triangle.H"
#include "oml/vector3d.h"
#include <iosfwd>

class RNLM
{
public:
    RNLM();
    // M&D eq 3.11 Alpha=Alpha_P, dR = CP^2 for nuclear
    // M&D eq 3.32 Alpha=Alpha_P*Alpha_Q/(Alpha_P+Alpha_Q) to electron repulsion
    //             dR = PQ^2
    RNLM(int Max, double Alpha, const RVec3& dR);
    virtual ~RNLM() {};
    bool CheckLMax(int L) const {return L<=itsLMax;}
    
    double operator()(const Polarization& p) const
    {
        return itsData(p.n,p.l,p.m);
    }

    // Element by element addition to this RNLM.
    void   Add  (const RNLM&, double Scale);
    void   Clear(                         )
    {
        itsData.Clear();
    }

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;
    virtual RNLM*    Clone(        ) const;

private:
    int          itsLMax;
    Triangle     itsData;
};

#endif //_RNLM_H_

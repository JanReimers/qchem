// File: RNLM.H  Manager the RNLM Auxillary functions.
#ifndef _RNLM_H_
#define _RNLM_H_



#include "Misc/Polarization.H"
#include "Misc/Triangle/Triangle.H"
#include "oml/vector3d.h"
#include "oml/imp/stream.h"

class RNLM
    : public virtual StreamableObject
{
public:
    RNLM();
    // M&D eq 3.11 Alpha=Alpha_P, dR = CP^2
    // M&D eq 3.32 Alpha=Alpha_P*Alpha_Q/(Alpha_P+Alpha_Q)
    //             dR = PQ^2
    RNLM(int Max, double Alpha, const RVec3& dR);

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
    Triangle     itsData;
};

#endif //_RNLM_H_

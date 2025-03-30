// File: Spin.H  A 3 state variable indicating the spin polarization.
#ifndef _Spin_H_
#define _Spin_H_

#include "Imp/Misc/pmstream.h"
#include <iosfwd>

class Spin : public virtual PMStreamableObject
{
public:
    enum State {Down,None,Up};

    Spin(       ) : itsState(None) {};
    Spin(State S) : itsState(S   ) {};

    bool IsPolarized  () const
    {
        return !(itsState==None);
    }
    int  GetDegeneracy() const
    {
        return IsPolarized() ? 1 : 2;
    }

    friend bool operator<(const Spin& a,const Spin& b)
    {
        int ia=static_cast<int>(a.itsState);
        int ib=static_cast<int>(b.itsState);
        return ia<ib;
    }
    friend bool operator==(const Spin& a,const Spin& b)
    {
        return a.itsState==b.itsState;
    }
    friend bool operator!=(const Spin& a,const Spin& b)
    {
        return a.itsState!=b.itsState;
    }

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&);
    size_t SequenceIndex() const //Used for op<
    {
        return static_cast<int>(itsState);
    }
    
    State itsState;
};

#endif //_Spin_H_

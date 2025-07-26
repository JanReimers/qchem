// File: BlochQN.H  A Quantum Number for atomic (spherical) symmetry
#ifndef _BlochQN_H_
#define _BlochQN_H_

//---------------------------------------------------------------------------------
//
//  Translational symmetry, Bloch function wave vector.
//

class BlochQN
    : public Symmetry
{
public:
    BlochQN(         );
    BlochQN(RVec3 theK);

    virtual int GetDegeneracy() const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual Symmetry* Clone(        ) const;

    RVec3   GetK() const
    {
        return itsK;
    }

private:
    RVec3 itsK;
};

#endif

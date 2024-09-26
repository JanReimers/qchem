// File: MeshImplementation.H   Mesh implemented as an vector of MeshPoint's.
#ifndef _MeshImplementation_H_
#define _MeshImplementation_H_



#include "Mesh/Mesh.H"
#include "oml/vector3d.h"
#include "oml/vector.h"
#include "oml/vector.h"

template <class T> class MeshIntegrator;
template <class T> class complex;

class MeshImplementation : public virtual Mesh
{
public:
    MeshImplementation(                 ) {};
    MeshImplementation(index_t numPoints);

    virtual void     ShiftOrigin (const RVec3&);
    virtual index_t  GetNumPoints() const
    {
        return itsPoints.size();
    }

    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

protected:
    void Initialize(const Vector<RVec3>& R,const Vector<double>& W);

private:
    friend class MeshBrowser;
    friend class MeshIntegrator<double>;
    friend class MeshIntegrator<std::complex<double> >;

    Vector<RVec3>  itsPoints;
    Vector<double> itsWeights;
};

#endif // _MeshImplementation_H_

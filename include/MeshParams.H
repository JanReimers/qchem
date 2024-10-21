// File: MeshParams.H Determine how integration meshes are generated.
#ifndef _MeshParams_H_
#define _MeshParams_H_

#include "types.H"

namespace qchem
{

enum RadialType {MHL,Log};
enum AngleType  {Gauss, GaussLegendre, EulerMclaren};


}

struct MeshParams
{
    qchem::RadialType radial_t;
    size_t     Nradial; //Number of mesh points in the radial directions
    size_t     MHL_m;       //m parameter for MHl radial mesh.
    double     MHL_alpha;
    // LoadFuzzyPoints also has an m paramater

    qchem::AngleType  angle_t;
    size_t     Nangle;    //Number of directionsGauss. 
    size_t     L_GL;;    //L for GaussLegendre, and EulerMclaren Nangle=(L+1)^2 /2
    size_t     m_GL;     //m parameter for GaussLegendre, and EulerMclaren angle meshes.
    
};


#endif // _MeshParams_H_

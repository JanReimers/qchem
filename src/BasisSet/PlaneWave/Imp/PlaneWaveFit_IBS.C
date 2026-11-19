// File: BasisSet/PlaneWave/Imp/PlaneWaveFit_IBS.C  The raster's voxel map under a crystal op (interface: ../PlaneWaveFit_IBS.C).
module;
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
module qchem.BasisSet.PlaneWave.PlaneWaveFit_IBS;

namespace qchem::BasisSet::PlaneWave
{

RasterVoxelMapT RasterVoxelMap(const Matrix3D<double>& W, const ivec3_t& N)
{
    const long n[3]={N.x, N.y, N.z};
    if (n[0]<=0 || n[1]<=0 || n[2]<=0) throw std::logic_error("RasterVoxelMap: an FFT grid needs N > 0 on every axis");
    RasterVoxelMapT M{};
    for (int a=0; a<3; ++a) for (int b=0; b<3; ++b)
    {
        const double m=W(a+1,b+1)*double(n[a])/double(n[b]);   // Matrix3D is 1-based
        M[a][b]=std::lround(m);
        if (std::fabs(m-double(M[a][b]))>1e-9)
            throw std::logic_error("RasterVoxelMap: the FFT grid "+std::to_string(n[0])+"x"+std::to_string(n[1])+"x"+std::to_string(n[2])
                                   +" is not invariant under a crystal op of the run's group (W_"+std::to_string(a+1)+std::to_string(b+1)
                                   +" N_"+std::to_string(a+1)+"/N_"+std::to_string(b+1)+" = "+std::to_string(m)+" is not an integer): "
                                   "the raster cannot be star-averaged on it -- the grid policy must hand an imposed run an invariant grid");
    }
    return M;
}

} //namespace qchem::BasisSet::PlaneWave

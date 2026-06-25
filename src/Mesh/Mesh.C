//  File: Mesh.C  Mesh integration PARAMETERS (the old AoS Mesh value type + MeshIntegrator are gone;
//  qcMesh1 owns the mesh and its quadrature now).  This module survives only to carry the typed
//  MeshParams that the DFT API (Ham_DFT*, CreateCD/VxcFitBasisSet, Fit_IBS::SetMesh) still passes.
module;
export module qchem.Mesh;
export import qchem.Types;

export namespace qchem
{
    //! \brief radial integration mesh types.
    enum RadialType {MHL,Log};
    //! \brief angular integration mesh types.
    enum AngleType  {Gauss, GaussLegendre, EulerMaclaren};
}

//! \brief parameters for controlling various Gaussian-quadrature integration mesh types.
//! Based on Murray, Handy & Laming, MOLECULAR PHYSICS 1993, 78, 997.
export struct MeshParams
{
    qchem::RadialType radial_t;
    size_t     Nradial;   //!<Number of mesh points in the radial directions
    size_t     MHL_m;     //!<m exponent for MHL radial mesh. m={1,2,3} are recommended.
    double     MHL_alpha; //!< \f$\alpha\f$ for the MHL radial mesh.  Related to atomic radius.
    qchem::AngleType  angle_t;
    size_t     Nangle;    //!<Number of directions for Gauss.
    size_t     L_GL;      //!<L for GaussLegendre, and EulerMclaren Nangle=(L+1)^2 /2
    size_t     m_GL;      //!<m parameter for GaussLegendre, and EulerMclaren angle meshes.
    //! Becke fuzzy-polyhedra smoothing parameter \f$m_{\mu}\f$ (4..12) for molecular meshes.
    size_t     m_mu;
};

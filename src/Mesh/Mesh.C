//  File: Mesh.C  mesh implementation
module;
#include <tuple>
#include <vector>

export module qchem.Mesh;
export import qchem.Types;

export class Mesh 
{ 
public:
    virtual ~Mesh () {};

    virtual void    ShiftOrigin(const RVec3&);
    virtual Mesh*   Clone      () const;
    
    typedef std::tuple<RVec3,double> rw_t;
    typedef std::vector<rw_t>         vec_t;
    typedef vec_t::const_iterator const_iterator;
    virtual const_iterator begin() const {return itsRWs.begin();}
    virtual const_iterator end  () const {return itsRWs.end  ();}
    virtual size_t         size () const {return itsRWs.size ();}

protected:
    void push_back(const RVec3& r, const double& w) {itsRWs.push_back(std::make_tuple(r,w));}
private:
    vec_t itsRWs;
};

export inline const RVec3 & r(const Mesh::rw_t& rw) {return std::get<0>(rw);}
export inline const double& w(const Mesh::rw_t& rw) {return std::get<1>(rw);}

export namespace qchem
{
    //! \brief radial integration mesh types.
    enum RadialType {MHL,Log};
    //! \brief angular integration mesh types.
    enum AngleType  {Gauss, GaussLegendre, EulerMaclaren};
}

//! \brief parameters for controlling various Gaussian Quadrature integration mesh types.
//! These are mostly based on the paper Quadrature schemes for integrals of density functional theory
//! Christopher W. Murray, Nicholas C. Handy & Gregory Laming, MOLECULAR PHVStCS, 1993, VOL. 78, No. 4, 997-1014\n 
//! Integrals are approximated as:\n
//! \f$\intop r^{2}dr\sin\theta d\theta d\phi F\left(r,\theta,\phi\right)\approx\sum_{i=1}^{N_{r}}w_{i}\sum_{j=1}^{N_{a}}w_{j}F\left(r_{i},\theta_{j},\phi_{j}\right)\f$ \n
//! where \f$ r_{i} \f$ and \f$ w_{i} \f$ are radial nodes and weights, and similarly for the angular integrals.
export struct MeshParams
{
    //! Supported radial mesh types are MHL and Log.  For MHL:\n
    //! \f$ r_{i}=\frac{\alpha\left(\frac{i}{N_{r}}\right)^{m}}{\left(1-\frac{i}{N_{r}}\right)^{m}} \f$ 
    //! \f$ w_{i}=\frac{\alpha mr_{i}^{2}}{N_{r}}\frac{\alpha\left(\frac{i}{N_{r}}\right)^{m-1}}{\left(1-\frac{i}{N_{r}}\right)^{m+1}} \f$ \n 
    //! And for the log mesh is defined by: \n
    //! \f$q_{r}=\frac{R_{max}}{R_{man}N_{r}},\quad q_{w}=\frac{1}{3}\left(\sqrt[3]{q_{r}}-\frac{1}{\sqrt[3]{q_{r}}}\right) \f$ \n 
    //! and \f$ r_{i}=q_{r}^{i},\quad w_{i}=r_{i}^{3}q_{w}\f$
    qchem::RadialType radial_t;
    size_t     Nradial; //!<Number of mesh points in the radial directions
    size_t     MHL_m;   //!<m exponent for MHl radial mesh. m={1,2,3} are recommended.
    double     MHL_alpha; //!< \f$\alpha\f$ for the MHL readial mesh.  Related to atomic radius.
    //! for angle_t=Gauss allowed values of Nangle are {1,2,6,8,12,24,30,32,50}.  For the GaussLegendre and
    //! EulerMclaren option spherical harmonics up to order L are integrated exactly.
    qchem::AngleType  angle_t;
    size_t     Nangle;    //!<Number of directions Gauss. 
    size_t     L_GL;     //!<L for GaussLegendre, and EulerMclaren Nangle=(L+1)^2 /2
    size_t     m_GL;     //!<m parameter for GaussLegendre, and EulerMclaren angle meshes.
    //! When combining multiple atom mesh basins into one molecular integration mesh one must
    //! blend the atom meshes with a weighting scheme as described by BECKE,A. D., 1988, J. chem. Phys., 88, 2547.
    //! In effect the Veroni polyhedra around each atom are continuously blended.  This is characterized
    //! another paramater \f$m_{\mu}\f$.  Recommeded ranges are \f$4\leq m_{\mu}\leq12\f$, with larger values
    //! being optimal for larger meshes.
    size_t     m_mu;
    
};



void Mesh::ShiftOrigin(const RVec3& r)
{
    for (auto& rw:itsRWs) std::get<0>(rw)+=r;
}

Mesh*   Mesh::Clone      () const
{
    return new Mesh(*this);
}

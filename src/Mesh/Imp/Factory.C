module;
#include <nlohmann/json.hpp>
module qchem.Mesh.Factory;
import qchem.Mesh.Internal.Types;
import qchem.Mesh.Internal.RadialTypes;
import oml;
using json = nlohmann::json;

namespace MeshF
{
    RadialMesh* Factory( qchem::RadialType rt,const nlohmann::json& js)
    {
        size_t N=js["N"];
        RadialMesh* m=0;
        switch (rt)
        {
        // case RadialType::Linear:
        //     {
        //         double start=js["start"].template get<double>(),stop=js["stop"].template get<double>();
        //         Vector3D<double> dir(0,0,1);//=js["dir"].template get<Vector3D<double>>();
        //         m=new LinearMesh(start,stop,dir,N);
        //     }
        //     break;
        case qchem::Log:
            {
                double start=js["start"].template get<double>(),stop=js["stop"].template get<double>();
                m=new LogRadialMesh(start,stop,N);
            }
            break;
        case qchem::MHL:
            {
                int mi=js["m"].template get<int>();
                double alpha=js["alpha"].template get<double>();
                m=new MHLRadialMesh(N,mi,alpha);
            }
            break;
        }
        assert(m);
        return m;
    }

    Mesh* Factory(qchem::AngleType at,const nlohmann::json& js)
    {
        Mesh* m=0;
        switch(at)
        {
            case qchem::EulerMaclaren:
            {
                int L=js["L"].template get<int>(),mi=js["m"].template get<int>();
                m=new EulerMaclarenAngularMesh(L,mi);
            }
            break;
            case qchem::Gauss:
            {
                int Nangle=js["Nangle"].template get<int>();
                m=new GaussAngularMesh(Nangle);
            }
            break;
            case qchem::GaussLegendre:
            {
                int L=js["L"].template get<int>(),mi=js["m"].template get<int>();
                m=new GaussLegendreAngularMesh(L,mi);
            }
            break;
        }
        assert(m);
        return m;
    }
}
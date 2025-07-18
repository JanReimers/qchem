module;
#include <nlohmann/json.hpp>
using json = nlohmann::json;

module qchem.Mesh.Factory;

import Mesh;
import Mesh.LinearMesh;
import Mesh.GaussAngularMesh;
import Mesh.GaussLegendreAngularMesh;
import Mesh.EulerMaclarenAngularMesh;
import Mesh.MHLRadialMesh;
import Mesh.LogRadialMesh;
import oml;

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
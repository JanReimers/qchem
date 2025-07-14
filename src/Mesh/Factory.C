// File: Mesh/Factory.C  Create various mesh types.

#include <Mesh/Factory.H>
#include "oml/vector3d.h"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

import Mesh;
import Mesh.LinearMesh;
import Mesh.GaussAngularMesh;
import Mesh.GaussLegendreAngularMesh;
import Mesh.EulerMaclarenAngularMesh;
import Mesh.MHLRadialMesh;
import Mesh.LogRadialMesh;

namespace MeshF
{
    RadialMesh* Factory( RadialType rt,const nlohmann::json& js)
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
        case RadialType::Log:
            {
                double start=js["start"].template get<double>(),stop=js["stop"].template get<double>();
                m=new LogRadialMesh(start,stop,N);
            }
            break;
        case RadialType::MHL:
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

    Mesh* Factory(AngularType at,const nlohmann::json& js)
    {
        Mesh* m=0;
        switch(at)
        {
            case AngularType::EulerMaclaren:
            {
                int L=js["L"].template get<int>(),mi=js["m"].template get<int>();
                m=new EulerMaclarenAngularMesh(L,mi);
            }
            break;
            case AngularType::Gauss:
            {
                int Nangle=js["Nangle"].template get<int>();
                m=new GaussAngularMesh(Nangle);
            }
            break;
            case AngularType::GaussLegendre:
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
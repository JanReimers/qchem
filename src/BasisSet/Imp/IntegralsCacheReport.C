// File: BasisSet/Imp/IntegralsCacheReport.C
module;
module qchem.BasisSet.IntegralsCacheReport;
import qchem.BasisSet.Internal.DB_Cache;   // theCache<T>() -- the mechanism, reachable here (same family) and nowhere public

namespace qchem::BasisSet
{
void EmitIntegralsCacheReport() { theCache<double>().EmitReport(); }
}

// File: Common/Imp/Parallel.C  The BLAS thread pin (rationale on the declaration).
module;
#include <cblas.h>   // openblas_set_num_threads -- an OpenBLAS extension, declared here by the
                     // openblas alternative of cblas.h (netlib's has no such call, which is the
                     // point: a BLAS swap breaks the build here rather than silently unpinning).
module qchem.Parallel;

namespace qchem {

void PinBlasToOneThread()
{
    openblas_set_num_threads(1);
}

} // namespace qchem

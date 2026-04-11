// File: LAParams.C Parameter for determing how the linear algebra for Lowden orthogonalization is done.
export module qchem.LAParams;

export namespace qchem 
{ 
    enum Ortho {Cholsky,Eigen,SVD}; //Basis set orthogonalization technique   
} //Namespace qchem

export struct LAParams
{
    qchem::Ortho  BasisOrthoAlgorithm;  //{Cholsky,Eigen,SVD}
    //For Ortho=SVD/Eigen we can truncate solutions with small
    //Eigen/singular values.  This removes linear dependency.
    double        TruncationTolerance;  
};

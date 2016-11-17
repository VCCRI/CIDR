#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::SelfAdjointEigenSolver;
using Rcpp::List;
using Rcpp::as;
using Rcpp::NumericMatrix;

// getEigenSpace - function to calculate eigenvalues & eigenvectors
// input: a NumericMatrix - cidr dissimilarity matrix
// returns: a lisf of values & vectors
// [[Rcpp::export]]
List getEigenSpace(NumericMatrix Xs) {
    const Map<MatrixXd>  X(as<Map<MatrixXd> >(Xs));
    SelfAdjointEigenSolver<MatrixXd> es(X);
    return List::create(Rcpp::Named("values") = es.eigenvalues(),
                        Rcpp::Named("vectors") = es.eigenvectors());
}

// eigen_centre - function to find product of input matrices
// - part of the calculation from the R package pcoa ; here we use the 
// c++ eigen library to calculate the product
// input: two NumericMatrix types - M & D
// returns: a list consisting of one entry, the result of the product M*D*M
// [[Rcpp::export]]
List eigen_centre(NumericMatrix Ms, NumericMatrix Ds) {
    const Map<MatrixXd>  M(as<Map<MatrixXd> >(Ms));
    const Map<MatrixXd>  D(as<Map<MatrixXd> >(Ds));
    
    return List::create(Rcpp::Named("centre") = M*D*M);
}

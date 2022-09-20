// Include guard
#ifndef __header_hpp__
#define __header_hpp__

#include<iomanip>
#include<cmath>
#include<armadillo>
#include<fstream>

// header file to be included in the main program project_2.cpp
// declaring functions defined in source.cpp

double max_offdiag_symm(const arma::mat& A, int& k, int& l);
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues,
  arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);
void max_offdiag_symm_test();

// End of include gurad
#endif

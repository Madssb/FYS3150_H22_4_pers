// Include guard
#ifndef __header_hpp__
#define __header_hpp__

#include<iomanip>
#include<cmath>
#include<armadillo>
#include<fstream>

using namespace std;

// header file to be included in the main program project_2.cpp
// declaring functions defined in source.cpp
double max_offdiag_symm(const arma::mat& A, int& k, int& l);

void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l);

void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues,
  arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged);

void max_offdiag_symm_test();

void ana_vs_arma_test(const double& a, const double& d, const double& e, const double& N);

arma::mat create_tridiag_mat(const double& a, const double& d, const double& e, const int& N);

void jacobi_eigensolver_multiple(const arma::mat& A, double eps, const int maxiter,
  int& iterations, bool& converged, fstream& outfile);

void three_lowest(const arma::vec& eigenvals, const arma::mat& eigenvecs);

void analytical_sol(const double& a, const double& d, const double& e, const double& N);

// End of include gurad
#endif

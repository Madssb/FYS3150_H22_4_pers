#include<iostream>
#include "header.hpp"

using namespace std;

int main(){

  // Testing Armadillo vs. analytical solution.
  int N = 6;
  double h = 1./6;
  double a = -1./(h*h);
  double d = 2./(h*h);
  double e = a;

  // ana_vs_arma_test(a, d, e, N);

  // Testing that the function max_offdiag_symm finds
  // the correct value.
  // max_offdiag_symm_test();

  // Testing the function jacobi_eigensolver and
  // jacobi_rotate.
  a = -2.5;
  d = 1.8;
  e = -2.5;
  arma::mat A = create_tridiag_mat(a, d, e, N);   // signature (-2,3,5)

  int iterations;
  arma::vec eigenvalues;
  arma::mat eigenvectors;

  double eps = 1e-8;
  int maxiter = 1e3;
  bool converged = false;

  jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter,
    iterations, converged);

  ana_vs_arma_test(a, d, e, N);

  return 0;
}

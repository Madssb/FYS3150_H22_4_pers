#include<iostream>
#include "header.hpp"

using namespace std;

int main(){

  // Testing Armadillo vs. analytical solution.
  int N_1 = 6;
  double h_1 = 1./7;
  double a_1 = -1./(h_1*h_1);
  double d_1 = 2./(h_1*h_1);
  double e_1 = a_1;

  string str_1 = "Testing Armadillo vs. analytical solution";
  string lines_1 = string(str_1.length(), '-');

  cout << endl << str_1 << endl;
  cout << lines_1 << endl;
  ana_vs_arma_test(a_1, d_1, e_1, N_1);
  cout << lines_1 << endl;

  // Testing that the function max_offdiag_symm finds
  // the correct value.
  string str_2 = "Testing that max_offiag_symm finds correct value";
  string lines_2 = string(str_2.length(), '-');

  cout << str_2 << endl;
  cout << lines_2 << endl;
  max_offdiag_symm_test();
  cout << lines_2 << endl;

  // Testing the function jacobi_eigensolver and
  // jacobi_rotate.
  double a_2 = -1/(h_1*h_1);
  double d_2 = 2/(h_1*h_1);
  double e_2 = -1/(h_1*h_1);
  arma::mat A_1 = create_tridiag_mat(a_2, d_2, e_2, N_1);

  int iterations = 0;
  arma::vec eigenvalues_1;
  arma::mat eigenvectors_1;

  double eps = 1e-8;
  int maxiter_1 = 1e5;
  bool converged = false;

  string str_3 = "Testing jacobi_eigensolver and jacobi_rotate";
  string lines_3 = string(str_3.length(), '-');

  cout << str_3 << endl;
  cout << lines_3 << endl;
  jacobi_eigensolver(A_1, eps, eigenvalues_1, eigenvectors_1, maxiter_1,
    iterations, converged);
  cout << lines_3 << endl;

  string str_4 = "Testing that the eigenvalues and eigenvectors from Armadillo\
  matches w/ analytical";
  string lines_4 = string(str_4.length(), '-');

  cout << str_4 << endl;
  cout << lines_4 << endl;
  ana_vs_arma_test(a_2, d_2, e_2, N_1);
  cout << lines_4 << endl;

  // Running the algorithm multiple times for different values of N, to
  // see how the number of iterations scales with N.
  int N_start = 2;
  int N_stop = 20;
  int maxiter_2 = 1e5;
  fstream outfile;

  // Comment out below to write to file.
  // for (int n = N_start; n <= N_stop; n++){
  //   bool converged = false;
  //   int iterations = 0;
  //
  //   arma::mat A_2 = create_tridiag_mat(a_2, d_2, e_2, n);
  //   jacobi_eigensolver_multiple(A_2, eps, maxiter_2, iterations, converged, outfile);
  // }

  // Solving the discrete equation for n = 10 and n = 100.
  for (int i = 0; i <= 1; i++){
    int n = 10 + (90*i);
    double h_2 = 1./n;
    int N_2 = n-1;

    double a_3 = -1./(h_2 * h_2);
    double d_3 = 2./(h_2 * h_2);
    double e_3 = -1./(h_2 * h_2);

    int iterations_3 = 0;
    bool converged_3 = false;

    arma::vec eigenvalues_2;
    arma::mat eigenvectors_2;

    arma::vec x_hat = arma::linspace(0, 1, n+1);
    arma::mat A_3 = create_tridiag_mat(a_3, d_3, e_3, N_2);

    jacobi_eigensolver(A_3, eps, eigenvalues_2, eigenvectors_2, maxiter_2,
    iterations_3, converged_3);

    // Uncomment below to write to file.
    // three_lowest(eigenvalues_2, eigenvectors_2);
  }
  return 0;
}

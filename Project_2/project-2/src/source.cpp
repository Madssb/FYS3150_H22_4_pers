#include "header.hpp"

using namespace std;

// Find the element w/ largest absolute value of matrix A.
// Since A is symmetric, we only scan through upper triangular
// part.
// Indices k and l are updated along the way.
double max_offdiag_symm(const arma::mat& A, int& k, int& l){

  int N = A.n_rows;
  double max = 0.;

  for (int i = 0; i < N; i++){

    for (int j = i+1; j < N; j++){

      if (abs(A(i,j)) > max){

        max = abs(A(i,j));
        k = i;
        l = j;
      }
    }
  }

  return max;
}

// Test function to test if max_offdiag_symm is working correctly.
void max_offdiag_symm_test(){

  arma::mat A = arma::mat(4, 4, arma::fill::eye);
  double max_offdiag;
  int k;
  int l;

  A(0, 3) = 0.5;
  A(1, 2) = -0.7;
  A(2, 1) = A(1, 2);
  A(3, 0) = A(0, 3);

  cout << endl << "Matrix A:" << endl;
  cout << A << endl;

  max_offdiag = max_offdiag_symm(A, k, l);

  cout << "Max off-diagonal absolute value:" << endl;
  cout << max_offdiag << endl;
  cout << "Element position:" << endl;
  cout << "A(" << k+1 << "," << l+1 << ")" << endl;
}

// Performs one rotation where the largest absolute value of
// A is rotated to zero.
void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){

  double tau = (A(l,l) - A(k,k)) / (2. * A(k,l));
  // Keeps track of elements so we dont use the wrong ones.
  double a_kk = A(k,k);
  double a_ll = A(l,l);
  double a_kl = A(k,l);
  double t;

  if (tau > 0){
    t = 1. / (tau + sqrt(1 + (tau * tau)));
  }
  else{
    t = -1. / (-tau + sqrt(1 + (tau * tau)));
  }

  double c = 1. / sqrt(1 + (t * t));
  double s = c * t;

  A(k,k) = a_kk - 2. * a_kl * c * s + a_ll * s * s;
  A(l,l) = a_ll + 2. * a_kl * c * s + a_kk * s * s;
  A(k,l) = 0.;
  A(l,k) = 0.;

  for (int i = 0; i < A.n_rows; i++){

    if (i != k && i != l){
      double a_ik = A(i,k);
      double a_il = A(i,l);

      A(i,k) = a_ik * c - a_il * s;
      A(k,i) = A(i,k);
      A(i,l) = a_il * c + a_ik * s;
      A(l,i) = A(i,l);
    }
  else{}
  }

  for (int i = 0; i < A.n_rows; i++){
    double r_ik = R(i,k);
    double r_il = R(i,l);

    R(i,k) = r_ik * c - r_il * s;
    R(i,l) = r_il * c + r_ik * s;
  }
}

// Runs jacobi_rotate until either
// 1. the iteration limit (maxiter) is reached
// 2. or the convergence criteria is reached.
void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues,
  arma::mat& eigenvectors, const int maxiter, int& iterations, bool& converged){

    int k;
    int l;

    arma::mat R = arma::mat(A.n_rows, A.n_cols);

    while (converged == false && iterations <= maxiter){
      double max_offdiagonal = max_offdiag_symm(A, k, l);

      iterations += 1;

      if (max_offdiagonal < eps){
        converged = true;
      }
      else{}
    }
    if (converged == false || iterations > maxiter){
      cout << endl << "Converged: " << converged << endl;
      cout << "Max iterations: " << maxiter << endl;
      cout << "Iterations: " << iterations << endl;
    }
    else{
      cout << endl << "Converged: " << converged << endl;

      eigenvalues = A.diag();
      eigenvectors = R.cols(arma::span::all);
    }
  }

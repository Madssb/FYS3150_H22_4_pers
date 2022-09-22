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

// Tests Armadillo solution and analytic solution of
// matrix eq. Ax = lambda x for tridiag. A = (NxN).
void ana_vs_arma_test(const double& a, const double& d, const double& e, const double& N){

  // Setting up tridiagonal matrix A w/ signature (a,d,e).
  arma::mat A = create_tridiag_mat(a, d, e, N);

  cout << endl << "Tridiagonal matrix A:" << endl;
  cout << A << endl;

  const double pi = 4. * atan(1.);

  arma::vec arma_eigenvals;
  arma::mat arma_eigenvecs;

  arma::vec ana_eigenvals = arma::vec(N);
  arma::mat ana_eigenvecs = arma::mat(N,N);

  for (int i = 0; i < N; i++){
    ana_eigenvals(i) = d + 2 * a * cos((i+1) * pi / (N+1));

    for (int n = 0; n < N; n++){
      double element = sin((n+1) * (i+1) * pi / (N+1));

      ana_eigenvecs.col(i)(n) = element;
    }
  }

  arma::eig_sym(arma_eigenvals, arma_eigenvecs, A);

  arma::vec ana_eigenvals_sorted = arma::sort(ana_eigenvals);
  arma::mat arma_eigenvecs_norm = arma::normalise(arma_eigenvecs, 1, 0);
  arma::mat ana_eigenvecs_norm = arma::normalise(ana_eigenvecs, 1, 0);

  cout << "EIGENVALUES" << endl;
  cout << "Armadillo:" << endl;
  cout << arma_eigenvals << endl;
  cout << "Analytical:" << endl;
  cout << ana_eigenvals_sorted << endl;

  cout << "EIGENVECTORS" << endl;
  cout << "Armadillo:" << endl;
  cout << arma_eigenvecs_norm << endl;
  cout << "Analytical:" << endl;
  cout << ana_eigenvecs_norm << endl;
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
  arma::mat R = arma::mat(4, 4, arma::fill::eye);
  int k_2=1;
  int l_2=2;
  jacobi_rotate(A, R, k_2, l_2);
  cout << A << endl;
  cout << R << endl;
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

  A(k,k) = a_kk * c * c - 2. * a_kl * c * s + a_ll * s * s;
  A(l,l) = a_ll * c * c + 2. * a_kl * c * s + a_kk * s * s;
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

    arma::mat R = arma::mat(A.n_rows, A.n_cols, arma::fill::eye);
    arma::mat A_m = A;

    while (converged == false && iterations <= maxiter){
      double max_offdiagonal = max_offdiag_symm(A_m, k, l);

      if (max_offdiagonal < eps){
        converged = true;
      }
      else{}

      jacobi_rotate(A_m, R, k, l);

      iterations += 1;
    }

    if (converged == false || iterations > maxiter){
      cout << endl << "Converged (1 is yes, 0 is no): " << converged << endl;
      cout << "Max iterations: " << maxiter << endl;
      cout << "Iterations: " << iterations << endl;
    }
    else{
      cout << endl << "Converged in " << iterations << " iterations for A=(" << A.n_rows <<
      "x" << A.n_rows << ")." << endl;

      eigenvalues = A_m.diag(0);
      eigenvectors = R.cols(arma::span::all);
    }
    if (A.n_rows <= 10){
      cout << endl << "Matrix A^(" << iterations << "):" << endl;
      // Printing out the matrix A_m where elements w/ value <= eps set to zero.
      cout << A_m.clean(eps) << endl;
    }
    else{
      cout << "Matrix is " << to_string(A.n_rows) << "x" << to_string(A.n_rows)
      << ", skipping print to terminal." << endl;
    }
  }

// Function to write N and no. of iterations to terminal
// for plotting in python
void jacobi_eigensolver_multiple(const arma::mat& A, double eps, const int maxiter,
  int& iterations, bool& converged, fstream& outfile){
    int k;
    int l;
    int N = A.n_rows;

    arma::mat R = arma::mat(N, N, arma::fill::eye);
    arma::mat A_m = A;

    while (converged == false && iterations <= maxiter){
      double max_offdiagonal = max_offdiag_symm(A_m, k, l);

      if (max_offdiagonal < eps){
        converged = true;
      }
      else{}

      jacobi_rotate(A_m, R, k, l);

      iterations += 1;
    }

    if (converged == false || iterations > maxiter){
      cout << endl << "Converged (1 is yes, 0 is no): " << converged << endl;
      cout << "Max iterations: " << maxiter << endl;
      cout << "Iterations: " << iterations << endl;
    }
    else{
      outfile.open("N_vs_iterations.txt", fstream::out | fstream::app);
      outfile << to_string(N) << "\t" << to_string(iterations) << endl;
      outfile.close();
    }
}

// Creating a tri-diagonal matrix
arma::mat create_tridiag_mat(const double& a, const double& d, const double& e, const int& N){

  arma::mat A = arma::mat(N,N);

  A.diag(-1).fill(a);
  A.diag(0).fill(d);
  A.diag(1).fill(e);

  return A;
}

// Sort 3 eigenvecs corresponding to the 3 lowest egenvals
// and write to .txt document
void three_lowest(const arma::vec& eigenvals, const arma::mat& eigenvecs){
  int N = eigenvals.size();

  arma::mat A = arma::mat(N+2, 4);
  arma::uvec indices = arma::sort_index(eigenvals);
  arma::vec x_hat = arma::linspace(0, 1, N+2);
  A.col(0) = x_hat;

  for (int i = 1; i < 4; i++){
    A(arma::span(1,N), i) = eigenvecs.col(indices(i));
  }
  string filename = "3_eigvecs_N_is_";
  string mat_size = to_string(N);
  string end_str = ".txt";
  string fullFilename = filename + mat_size + end_str;
  fstream outfile;

  outfile.open(fullFilename, fstream::out | fstream::app);
  outfile << A << endl;
  outfile.close();
}

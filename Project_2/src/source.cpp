#include "header.hpp"
#include <iostream>
using namespace std;

// Find the element w/ largest absolute value of matrix A.
// Since A is symmetric, we only scan through upper triangular
// part.
// Indices k and l are updated along the way.
double max_offdiag_symm(const arma::mat &symmetric_matrix, int &k, int &l)
{

  int N = symmetric_matrix.n_rows;
  double max = 0.;

  for (int i = 0; i < N; i++)
  {

    for (int j = i + 1; j < N; j++)
    {

      if (abs(symmetric_matrix(i, j)) > max)
      {

        max = abs(symmetric_matrix(i, j));
        k = i;
        l = j;
      }
    }
  }

  return max;
}

// Tests Armadillo solution and analytic solution of
// matrix eq. Ax = lambda x for tridiag. A = (NxN).
void ana_vs_arma_test(const double &a, const double &d, const double &e, const double &N)
{

  // Setting up tridiagonal matrix A w/ signature (a,d,e).
  arma::mat A = create_tridiag_mat(a, d, e, N);

  cout << std::endl
       << "Tridiagonal matrix A:" << std::endl;
  cout << A << std::endl;

  const double pi = 4. * atan(1.);

  arma::vec arma_eigenvals;
  arma::mat arma_eigenvecs;

  arma::vec ana_eigenvals = arma::vec(N);
  arma::mat ana_eigenvecs = arma::mat(N, N);

  for (int i = 0; i < N; i++)
  {
    ana_eigenvals(i) = d + 2 * a * cos((i + 1) * pi / (N + 1));

    for (int n = 0; n < N; n++)
    {
      double element = sin((n + 1) * (i + 1) * pi / (N + 1));

      ana_eigenvecs.col(i)(n) = element;
    }
  }

  arma::eig_sym(arma_eigenvals, arma_eigenvecs, A);

  arma::vec ana_eigenvals_sorted = arma::sort(ana_eigenvals);
  arma::mat arma_eigenvecs_norm = arma::normalise(arma_eigenvecs, 1, 0);
  arma::mat ana_eigenvecs_norm = arma::normalise(ana_eigenvecs, 1, 0);

  cout << "EIGENVALUES" << std::endl;
  cout << "Armadillo:" << std::endl;
  cout << arma_eigenvals << std::endl;
  cout << "Analytical:" << std::endl;
  cout << ana_eigenvals_sorted << std::endl;

  cout << "EIGENVECTORS" << std::endl;
  cout << "Armadillo:" << std::endl;
  cout << arma_eigenvecs_norm << std::endl;
  cout << "Analytical:" << std::endl;
  cout << ana_eigenvecs_norm << std::endl;

  cout << "RELATIVE ERRORS" << std::endl;
  cout << "Eigenvalues:" << std::endl;
  cout << abs(ana_eigenvals_sorted - arma_eigenvals) / abs(arma_eigenvals) << std::endl;
  cout << "Eigenvectors:" << std::endl;
  cout << abs(ana_eigenvecs_norm - arma_eigenvecs_norm) / abs(arma_eigenvecs_norm) << std::endl;
}

// Test function to test if max_offdiag_symm is working correctly.
void max_offdiag_symm_test()
{

  arma::mat A = arma::mat(4, 4, arma::fill::eye);
  double max_offdiag;
  int k;
  int l;

  A(0, 3) = 0.5;
  A(1, 2) = -0.7;
  A(2, 1) = A(1, 2);
  A(3, 0) = A(0, 3);

  cout << std::endl
       << "Matrix A:" << std::endl;
  cout << A << std::endl;

  max_offdiag = max_offdiag_symm(A, k, l);

  cout << "Max off-diagonal absolute value:" << std::endl;
  cout << max_offdiag << std::endl;
  cout << "Element position:" << std::endl;
  cout << "A(" << k + 1 << "," << l + 1 << ")" << std::endl;
  arma::mat R = arma::mat(4, 4, arma::fill::eye);
  int k_2 = 1;
  int l_2 = 2;
  jacobi_rotate(A, R, k_2, l_2);
}

// Performs one rotation where the largest absolute value of
// A is rotated to zero.
void jacobi_rotate(arma::mat &symmetric_matrix, arma::mat &R, int k, int l)
{

  double tau = (symmetric_matrix(l, l) - symmetric_matrix(k, k)) / (2. * symmetric_matrix(k, l));
  // Keeps track of elements so we dont use the wrong ones.
  double a_kk = symmetric_matrix(k, k);
  double a_ll = symmetric_matrix(l, l);
  double a_kl = symmetric_matrix(k, l);
  double t;

  if (tau > 0)
  {
    t = 1. / (tau + sqrt(1 + (tau * tau)));
  }
  else
  {
    t = -1. / (-tau + sqrt(1 + (tau * tau)));
  }

  double c = 1. / sqrt(1 + (t * t));
  double s = c * t;

  symmetric_matrix(k, k) = a_kk * c * c - 2. * a_kl * c * s + a_ll * s * s;
  symmetric_matrix(l, l) = a_ll * c * c + 2. * a_kl * c * s + a_kk * s * s;
  symmetric_matrix(k, l) = 0.;
  symmetric_matrix(l, k) = 0.;

  for (int i = 0; i < symmetric_matrix.n_rows; i++)
  {

    if (i != k && i != l)
    {
      double a_ik = symmetric_matrix(i, k);
      double a_il = symmetric_matrix(i, l);

      symmetric_matrix(i, k) = a_ik * c - a_il * s;
      symmetric_matrix(k, i) = symmetric_matrix(i, k);
      symmetric_matrix(i, l) = a_il * c + a_ik * s;
      symmetric_matrix(l, i) = symmetric_matrix(i, l);
    }
    else
    {
    }
  }

  for (int i = 0; i < symmetric_matrix.n_rows; i++)
  {
    double r_ik = R(i, k);
    double r_il = R(i, l);

    R(i, k) = r_ik * c - r_il * s;
    R(i, l) = r_il * c + r_ik * s;
  }
}

// Runs jacobi_rotate until either
// 1. the iteration limit (maxiter) is reached
// 2. or the convergence criteria is reached.
void jacobi_eigensolver(const arma::mat &A, double eps, arma::vec &eigenvalues,
                        arma::mat &eigenvectors, const int maxiter, int &iterations, bool &converged)
{

  int k;
  int l;

  arma::mat R = arma::mat(A.n_rows, A.n_cols, arma::fill::eye);
  arma::mat A_m = A;

  while (converged == false && iterations <= maxiter)
  {
    double max_offdiagonal = max_offdiag_symm(A_m, k, l);

    if (max_offdiagonal < eps)
    {
      converged = true;
    }
    else
    {
    }

    jacobi_rotate(A_m, R, k, l);

    iterations += 1;
  }

  if (converged == false || iterations > maxiter)
  {
    cout << std::endl
         << "Converged (1 is yes, 0 is no): " << converged << std::endl;
    cout << "Max iterations: " << maxiter << std::endl;
    cout << "Iterations: " << iterations << std::endl;
  }
  else
  {
    cout << std::endl
         << "Converged in " << iterations << " iterations for A=(" << A.n_rows << "x" << A.n_rows << ")." << std::endl;

    eigenvalues = A_m.diag(0);
    eigenvectors = R.cols(arma::span::all);
  }
  if (A.n_rows <= 10)
  {
    cout << std::endl
         << "Matrix A^(" << iterations << "):" << std::endl;
    // Printing out the matrix A_m where elements w/ value <= eps set to zero.
    cout << A_m.clean(eps) << std::endl;
  }
  else
  {
    cout << "Matrix is " << to_string(A.n_rows) << "x" << to_string(A.n_rows)
         << ", skipping print to terminal." << std::endl;
  }
}

//
// Function to write N and # of iterations to terminal
// for plotting in python
//

void jacobi_eigensolver_multiple(const arma::mat &A, double eps, const int maxiter,
                                 int &iterations, bool &converged, fstream &outfile)
{
  int k;
  int l;
  int N = A.n_rows;

  arma::mat R = arma::mat(N, N, arma::fill::eye);
  arma::mat A_m = A;

  while (converged == false && iterations <= maxiter)
  {
    double max_offdiagonal = max_offdiag_symm(A_m, k, l);

    if (max_offdiagonal < eps)
    {
      converged = true;
    }
    else
    {
    }

    jacobi_rotate(A_m, R, k, l);

    iterations += 1;
  }

  if (converged == false || iterations > maxiter)
  {
    cout << std::endl
         << "Converged (1 is yes, 0 is no): " << converged << std::endl;
    cout << "Max iterations: " << maxiter << std::endl;
    cout << "Iterations: " << iterations << std::endl;
  }
  else
  {
    outfile.open("N_vs_iterations.txt", fstream::out | fstream::app);
    outfile << to_string(N) << "\t" << to_string(iterations) << std::endl;
    outfile.close();
  }
}

// Creating a tri-diagonal matrix
arma::mat create_tridiag_mat(const double &a, const double &d, const double &e, const int &N)
{

  arma::mat A = arma::mat(N, N);

  A.diag(-1).fill(a);
  A.diag(0).fill(d);
  A.diag(1).fill(e);

  return A;
}

// Sort 3 eigenvecs corresponding to the 3 lowest egenvals
// and write to .txt document
void three_lowest(const arma::vec &eigenvals, const arma::mat &eigenvecs)
{
  int N = eigenvals.size();

  arma::mat A = arma::mat(N + 2, 4);
  arma::uvec indices = arma::sort_index(eigenvals);
  arma::vec x_hat = arma::linspace(0, 1, N + 2);
  A.col(0) = x_hat;

  for (int i = 1; i < 4; i++)
  {
    A(arma::span(1, N), i) = eigenvecs.col(indices(i));
  }
  string filename = "3_eigvecs_N_is_";
  string mat_size = to_string(N);
  string end_str = ".txt";
  string fullFilename = filename + mat_size + end_str;
  fstream outfile;

  outfile.open(fullFilename, fstream::out | fstream::app);
  outfile << A << std::endl;
  outfile.close();
}

// Analytical eigenvalues and eigenvectors
void analytical_sol(const double &a, const double &d, const double &e, const double &N)
{
  arma::mat A = create_tridiag_mat(a, d, e, N);
  arma::mat B = arma::mat(N + 2, 3);
  const double pi = 4. * atan(1.);

  arma::vec eigenvals = arma::vec(N);
  arma::mat eigenvecs = arma::mat(N, N);

  for (int i = 0; i < N; i++)
  {
    eigenvals(i) = d + 2 * a * cos((i + 1) * pi / (N + 1));

    for (int n = 0; n < N; n++)
    {
      double element = sin((n + 1) * (i + 1) * pi / (N + 1));

      eigenvecs.col(i)(n) = element;
    }
  }
  arma::uvec indices = arma::sort_index(eigenvals);

  for (int i = 0; i < 3; i++)
  {
    B(arma::span(1, N), i) = eigenvecs.col(indices(i));
  }
  string filename = "3_ANA_eigvecs_N_is_";
  string mat_size = to_string(N);
  string end_str = ".txt";
  string fullFilename = filename + mat_size + end_str;
  fstream outfile;
  outfile.open(fullFilename, fstream::out | fstream::app);
  outfile << B << std::endl;
  outfile.close();
}

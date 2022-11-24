/*
This file contains the definition of the Matrix class
*/

#include "matrix.hpp"
#include <complex>

// Constructor that takes no. of steps in space (M) and time (N), and the
// total time (T)
Matrix::Matrix(double h_in, double dt_in, double T_in)
{
  h = h_in;
  dt = dt_in;
  T = T_in;

  M = 1. / h + 1;
  N = 1. / dt + 1;
  r = arma::cx_double(0, dt / (2 * h * h));

  U.zeros(M - 2, M - 2);
  V.zeros(M - 2, M - 2);
  S.zeros(M - 2, M - 2, N);

  A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
  B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

// Get index k from (i, j)
int Matrix::index_k(int i, int j)
{
  return i + (M - 2) * j;
}

void Matrix::set_potential(std::string filename, double V_0)
{

}

void Matrix::fill_matrices()
{
  // Vectors a and b, which will be the main diagonals of A and B
  arma::cx_vec a((M - 2) * (M - 2), arma::fill::zeros);
  arma::cx_vec b((M - 2) * (M - 2), arma::fill::zeros);

  // Tridiagonal sub matrix with signature (r, 0, r)
  arma::cx_mat sub_diag(M - 2, M - 2, arma::fill::zeros);
  sub_diag.diag(-1).fill(r);
  sub_diag.diag(1).fill(r);

  // Filling diagonals of A and B
  A.diag(M - 2).fill(-r);
  A.diag(2 - M).fill(-r);
  B.diag(M - 2).fill(r);
  B.diag(2 - M).fill(r);

  // Filling a and b
  arma::cx_double a_k, b_k;
  int k;

  for (int i = 0; i < M - 2; i++)
  {
    for (int j = 0; j < M - 2; j++)
    {
      k = index_k(i, j);

      a_k = arma::cx_double(1., 4. * r.imag() + .5 * dt * V(i, j));
      b_k = std::conj(a_k);

      a(k) = a_k;
      b(k) = b_k;
    }
  }

  // Filling A and B with sub matrices
  for (int i = 0; i < M - 2; i++)
  {
    int j = i * (M - 2);

    A.submat(j, j, j + M - 3, j + M - 3) = -sub_diag;
    B.submat(j, j, j + M - 3, j + M - 3) = sub_diag;
  }

  // Setting main diagonals of A and B
  A.diag() = a;
  B.diag() = b;
}

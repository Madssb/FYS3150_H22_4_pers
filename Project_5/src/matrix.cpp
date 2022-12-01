/*
This file contains the definition of the Matrix class
*/

#include "matrix.hpp"
#include <complex>
#include <fstream>
#include <assert.h>

// Constructor takes spatial and time step lentgth as well as total time
Matrix::Matrix(double h_in, double dt_in, double T_in)
{
  h = h_in;
  dt = dt_in;
  T = T_in;

  M = 1. / h + 1;
  N = T / dt + 1;
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

// Set potential barriers
void Matrix::set_potential(std::string filename)
{
  // Load the file into the potential
  V.load(filename, arma::raw_ascii);
  assert(V.n_rows == (M - 2) && V.n_cols == (M - 2));

}

// Fill the matrices according to the Crank-Nicolson regime
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

// Solve matrix eq. Au^(n+1) = Bu^n
void Matrix::solve()
{
  // Solving the equation for each time step and saving as slice in S
  for (int i = 1; i < N; i++)
  {
    u_current = U.as_col();
    u_new = arma::spsolve(A, B * u_current);

    // Refilling u vector into U matrix
    for (int n = 0; n < U.n_cols; n ++)
    {
      for (int k = 0; k < M - 2; k++)
      {
        U.col(n)(k) = u_new(k + n * (M - 2));
      }
    }

    S.slice(i) = U;
    u_current = u_new;
  }
}

// Set the initial state of the system
void Matrix::set_initial_state(double x_c, double sigma_x, double p_x, double y_c,
                       double sigma_y, double p_y)
{
  arma::cx_double exponent;
  double re, im;
  int xi, yi;

  for (int i = 0; i < M - 2; i++)
  {
    yi = i;

    for (int j = 0; j < M - 2; j++)
    {
      xi = j;

      re = -(xi * h - x_c) * (xi * h * x_c) / (2 * sigma_x * sigma_x)
           -(yi * h - y_c) * (yi * h - y_c) / (2 * sigma_y * sigma_y);

      im = p_x * (xi * h - x_c) + (yi * h - y_c);

      exponent = arma::cx_double(re, im);

      U(i, j) = std::exp(exponent);
    }
  }

  // Normalize
  U /= std::sqrt(arma::accu(arma::conj(U) % U));
  S.slice(0) = U;
}

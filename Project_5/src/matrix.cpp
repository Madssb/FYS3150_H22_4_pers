/*
This file contains the definition of the Matrix class
*/

#include "matrix.hpp"
#include <complex>
#include <fstream>

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
  // S.zeros(M - 2, M - 2, N);

  u_new = U.as_col();
  u_current = U.as_col();

  A.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
  B.zeros((M - 2) * (M - 2), (M - 2) * (M - 2));
}

// Get index k from (i, j)
int Matrix::index_k(int i, int j)
{
  return i + (M - 2) * j;
}

// Set potential barriers
void Matrix::set_potential(std::string filename, double V_0)
{
  // Extracting potenital parameters from file
  arma::rowvec params;
  params.load(filename, arma::raw_ascii);

  int n_slits = params(0);        // No.  of slits
  double x_pos = params(1);       // Position (center) at the x axis
  double x_width = params(2);     // Width of barrier along x axis
  double app = params(3);         // Apperture
  double c_length = params(4);    // Length of (one of the) center barrier(s)

  int n_center = n_slits - 1; // No. of center barriers
  double e_length = .5 * (1. - c_length - n_slits * app);

  double xi, yi;

  for (int i = 0; i < M - 2; i++)
  {
    yi = (i + 1) * h;

    for (int j = 0; j < M - 2; j++)
    {
      xi = (j + 1) * h;

      if (xi >= x_pos - .5 * x_width && xi <= x_pos + .501 * x_width)
      {
        if (yi <= e_length)
        {
          V(i, j) = V_0;
        }
        else if (yi >= e_length + app && yi <= e_length + app + c_length)
        {
          V(i, j) = V_0;
        }
        else if (yi >= 1. - e_length)
        {
          V(i, j) = V_0;
        }
      }
    }
  }
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
void Matrix::solve(std::string out_filename)
{
  std::ofstream outfile;
  outfile.open(out_filename);

  for (int i = 0; i < N; i++)
  {
    outfile << U << std::endl;
    // S.slice(i) = U;
    u_new = arma::spsolve(A, B * u_current);
    u_current = u_new;
  }
  
  outfile.close();
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

      re = ((xi - x_c) * (xi - x_c)) / (2 * sigma_x * sigma_x)
         - ((yi - y_c) * (yi - y_c)) / (2 * sigma_y * sigma_y);

      im = p_x * (xi - x_c) + (yi - y_c);

      exponent = arma::cx_double(re, im);

      U(i, j) = exp(exponent);
    }
  }
}

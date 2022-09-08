// Task 8
#include <iostream>
#include <armadillo>
#include <cmath>
#include <iomanip>

// no. of points to simulate
int n = 10000;

// step length
double h = 1./n;

// elements of diagonals
double a = -1.;
double b = 2.;
double c = -1.;

// source term
double f(double x)
{
  double res = 100. * std::exp(-10. * x);
  return res;
}

// analytic solution
double u_sol(double x)
{
  double u = 1. - (1. - std::exp(-10.)) * x - std::exp(-10. * x);
  return u;
}

// g solution
double g_sol(double x)
{
  double g = h * h * (f(x) - a/b * f(x));
  return g;
}

// declare vectors
arma::vec x = arma::linspace(0, 1, n+1);
arma::vec u = arma::vec(n+1);
arma::vec v = arma::vec(n+1);
arma::vec b_tilde = arma::vec(n+1);
arma::vec g_tilde = arma::vec(n+1);

int main()
{
  // initial values
  double g_0 = g_sol(x(0));
  double u_0 = u_sol(x(0));
  double b_0 = b;

  // declare variables to be used in loop
  double g_i;
  double w_i;

  // initialize vectors
  u(0) = u_0;
  b_tilde(0) = b_0;
  g_tilde(0) = g_0;

  // forward sweep
  for (int i = 1; i <= n; i++)
  {
    g_i = g_sol(x(i));
    w_i = a/b_tilde(i-1);

    u(i) = u_sol(x(i));

    b_tilde(i) = b - w_i * c;
    g_tilde(i) = g_i - w_i * g_tilde(i-1);
  }

  v(n) = g_tilde(n) / b_tilde(n);

  // backwards sweep
  for (int i = n; i > 0; i--)
  {
    v(i-1) = (g_tilde(i-1) - c * v(i) / b_tilde(i-1));
  }

  x.save("x.bin");
  u.save("u.bin");
  v.save("v.bin");

  return 0;
}

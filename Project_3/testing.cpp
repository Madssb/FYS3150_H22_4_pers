// Code for testing functionality
#include "particle.hpp"
#include "penning_trap.hpp"
#include <cmath>
#include <complex>

using namespace std;

/// Defining constants
const double T = 9.64852558 * 10;   // [u/(mu*s)e ]
const double V = 9.64852558 * 1e7;  // [u(mu*m)^2/(mu*s)^2e]
const double B_0 = 1. * T;          // Magnetic field constant
const double V_0 = 10. * V;         // Electric field constant
const double d = 1e4;               // [mu*m]

// Analytical solution
void analytical_solution(double q, double m, double x_0, double v_0,
                              double z_0, arma::vec t, arma::mat& r)
{
  const double w_0 = q * B_0 / m;
  double w_z = 2 * q * V_0 / (m * pow(d, 2.));

  const double w_plus = (w_0 + sqrt(pow(w_0, 2.) - 2 * w_z)) / 2.;
  const double w_min = (w_0 - sqrt(pow(w_0, 2.) - 2 * w_z)) / 2.;

  const double A_plus = (v_0 + w_min * x_0) / (w_min - w_plus);
  const double A_min = - (v_0 + w_min * x_0) / (w_min - w_plus);
  const complex<double> i(0,1);

  r.row(0) = arma::vec{x_0, 0., z_0}.t();

  int N = t.size();

  for (int i = 1; i < N; i++)
  {
    double z = z_0 * cos(sqrt(w_z) * t(i));
    complex<double> f = A_plus * exp(- i * w_plus * t(i))
                      + A_min * exp(- i * w_min * t(i));

    double x = f.real();
    double y = f.imag();

    r.row(i) = arma::vec{x, y, z}.t();
  }
}


int main()
{

  PenningTrap trap = PenningTrap(B_0, V_0, d);

  double q = 1.;
  double m = 1.;

  double x_0 = 1, y_0 = 1, z_0 = 1;
  double vx_0 = 0, vy_0 = 0, vz_0 = 0;
  double v_0 = 0;

  arma::vec r_0 = arma::vec{x_0, y_0, z_0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r_0, v);

  trap.add_particle(p);

  int N = 10;
  double dt = 0.1;
  double steps = N / dt;

  arma::vec t = arma::linspace(0, N, steps);
  arma::mat r = arma::mat(steps, 3);
  arma::mat sol = arma::mat(steps, 4);

  analytical_solution(q, m, x_0, v_0, z_0, t, r);

  sol.col(0) = t;
  sol.cols(1, 3) = r;

  cout << sol << endl;

  return 0;
}

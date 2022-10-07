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
double analytical_solution_z(double q, double V_0, double m, double d, double t)
{
  double w_z_squared = 2 * q * V_0 / (m * pow(d, 2.));
  double w_z = sqrt(w_z_squared);
  double z_0 = 1.;
  double z = z_0 * cos(w_z * t);

  return z;
}

int main()
{
  int width = 18;
  int prec = 4;

  PenningTrap trap = PenningTrap(B_0, V_0, d);

  double q = 1.;
  double m = 1.;

  double x_0 = 1, y_0 = 0, z_0 = 1;
  double vx_0 = 0, vy_0 = 1, vz_0 = 0;

  arma::vec r_0 = arma::vec{x_0, y_0, z_0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r_0, v);

  trap.add_particle(p);

  int N = 10;
  double dt = 1e-3;
  int steps = N / dt;

  arma::vec t = arma::linspace(0, N, steps);

  cout << "#" << setw(width) << "Time"
       << setw(width) << "Num"
       << setw(width) << "Ana"
       << endl;

  for (int i = 0; i < steps; i++)
  {
    trap.evolve_RK4(dt);

    double z_num = trap.particles[0].r[2];
    double z_ana = analytical_solution_z(p.q, V_0, p.m, d, t(i));

    cout << setw(width) << setprecision(prec) << t(i)
         << setw(width) << setprecision(prec) << z_num
         << setw(width) << setprecision(prec) << z_ana
         << endl;
  }

  return 0;
}

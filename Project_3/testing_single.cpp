// Code for testing functionality
#include "particle.hpp"
#include "penning_trap.hpp"
#include "analytical_z.hpp"
#include <cmath>

using namespace std;

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

  int N = 100;
  double dt = 1e-3;
  int steps = N / dt;

  arma::vec t = arma::linspace(0, N, steps);

  cout << "#" << setw(width) << "Time"
       << setw(width) << "Num"
       << setw(width) << "Ana"
       << endl;

  for (int i = 0; i < steps; i++)
  {
    trap.evolve_RK4(dt, false);

    double z_num = trap.particles[0].r[2];
    double z_ana = analytical_solution_z(p.q, V_0, p.m, d, t(i));

    cout << setw(width) << setprecision(prec) << t(i)
         << setw(width) << setprecision(prec) << z_num
         << setw(width) << setprecision(prec) << z_ana
         << endl;
  }

  return 0;
}

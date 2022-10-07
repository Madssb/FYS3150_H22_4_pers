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

  arma::vec r_1 = arma::vec{x_0, y_0, z_0};
  arma::vec r_2 = - r_1;
  arma::vec v_1 = arma::vec{vx_0, vy_0, vz_0};
  arma::vec v_2 = 2. * v_1;

  Particle p_1 = Particle(q, m, r_1, v_1);
  Particle p_2 = Particle(q, m, r_2, v_2);

  trap.add_particle(p_1);
  trap.add_particle(p_2);

  int N = 3;
  double dt = 1e-4;
  int steps = N / dt;

  arma::vec t = arma::linspace(0, N, steps);

  cout << "#" << setw(width) << "Time"
       << setw(width) << "x1"
       << setw(width) << "x2"
       << setw(width) << "y1"
       << setw(width) << "y2"
       << endl;

  for (int i = 0; i < steps; i++)
  {
    trap.evolve_RK4(dt, true);

    double x1 = trap.particles[0].r[0];
    double x2 = trap.particles[1].r[0];
    double y1 = trap.particles[0].r[1];
    double y2 = trap.particles[1].r[1];


    cout << setw(width) << setprecision(prec) << t(i)
         << setw(width) << setprecision(prec) << x1
         << setw(width) << setprecision(prec) << x2
         << setw(width) << setprecision(prec) << y1
         << setw(width) << setprecision(prec) << y2
         << endl;
  }

  return 0;
}

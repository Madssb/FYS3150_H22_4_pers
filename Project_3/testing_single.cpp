// Code for testing functionality
#include "particle.hpp"
#include "penning_trap.hpp"
#include "analytical.hpp"
#include <cmath>

using namespace std;

int main()
{
  int width = 18;
  int prec = 4;

  PenningTrap trap = PenningTrap(B_0, V_0, d);

  double q = 1.;
  double m = 40.078;  // Ca+ ion

  double x_0 = 1., y_0 = 0., z_0 = 1.;
  double vx_0 = 0., vy_0 = 1., vz_0 = 0.;

  arma::vec r = arma::vec{x_0, y_0, z_0};
  arma::vec v = arma::vec{vx_0, vy_0, vz_0};

  Particle p = Particle(q, m, r, v);

  trap.add_particle(p);

  int total_time = 100;
  double dt = 1e-3;
  int steps = total_time / dt;

  arma::vec t = arma::linspace(0, total_time, steps);

  cout << "#" << setw(width) << "Time"
       << setw(width) << "x_num"
       << setw(width) << "y_num"
       << setw(width) << "z_num"
       << setw(width) << "x_ana"
       << setw(width) << "y_ana"
       << setw(width) << "z_ana"
       << endl;

  cout << setw(width) << setprecision(prec) << t(0)
       << setw(width) << setprecision(prec) << x_0
       << setw(width) << setprecision(prec) << y_0
       << setw(width) << setprecision(prec) << z_0
       << setw(width) << setprecision(prec) << x_0
       << setw(width) << setprecision(prec) << y_0
       << setw(width) << setprecision(prec) << z_0
       << endl;


  for (int i = 1; i < steps; i++)
  {
    trap.evolve_RK4(dt, true);

    arma::vec r_ana = analytical_solution(r, v, p.q, B_0, p.m, V_0, d, t(i));

    double x_ana = r_ana[0], y_ana = r_ana[1], z_ana = r_ana[2];
    double x_num = trap.particles[0].r[0];
    double y_num = trap.particles[0].r[1];
    double z_num = trap.particles[0].r[2];

    cout << setw(width) << setprecision(prec) << t(i)
         << setw(width) << setprecision(prec) << x_num
         << setw(width) << setprecision(prec) << y_num
         << setw(width) << setprecision(prec) << z_num
         << setw(width) << setprecision(prec) << x_ana
         << setw(width) << setprecision(prec) << y_ana
         << setw(width) << setprecision(prec) << z_ana
         << endl;
  }

  return 0;
}

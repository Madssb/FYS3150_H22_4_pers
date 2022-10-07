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

  double x_0 = 1., y_0 = 0., z_0 = 1.;
  double vx_0 = 0., vy_0 = 1., vz_0 = 0.;

  arma::vec r_1 = arma::vec{x_0, y_0, z_0};
  arma::vec r_2 = - r_1;
  arma::vec v_1 = arma::vec{vx_0, vy_0, vz_0};
  arma::vec v_2 = 2. * v_1;

  Particle p_1 = Particle(q, m, r_1, v_1);
  Particle p_2 = Particle(q, m, r_2, v_2);

  trap.add_particle(p_1);
  trap.add_particle(p_2);

  int total_time = 100;   // [mu*s]
  double dt = 1e-3;
  int steps = total_time / dt;

  arma::vec t = arma::linspace(0, total_time, steps);

  cout << "#" << setw(width) << "Time"
       << setw(width) << "x1"
       << setw(width) << "x2"
       << setw(width) << "y1"
       << setw(width) << "y2"
       << setw(width) << "z1"
       << setw(width) << "z2"
       << setw(width) << "vx1"
       << setw(width) << "vx2"
       << setw(width) << "vy1"
       << setw(width) << "vy2"
       << setw(width) << "vz1"
       << setw(width) << "vz2"
       << endl;

  for (int i = 0; i < steps; i++)
  {
    trap.evolve_RK4(dt, true);

    double x1 = trap.particles[0].r[0];
    double x2 = trap.particles[1].r[0];
    double y1 = trap.particles[0].r[1];
    double y2 = trap.particles[1].r[1];
    double z1 = trap.particles[0].r[2];
    double z2 = trap.particles[1].r[2];
    double vx1 = trap.particles[0].v[0];
    double vx2 = trap.particles[1].v[0];
    double vy1 = trap.particles[0].v[1];
    double vy2 = trap.particles[1].v[1];
    double vz1 = trap.particles[0].v[2];
    double vz2 = trap.particles[1].v[2];


    cout << setw(width) << setprecision(prec) << t(i)
         << setw(width) << setprecision(prec) << x1
         << setw(width) << setprecision(prec) << x2
         << setw(width) << setprecision(prec) << y1
         << setw(width) << setprecision(prec) << y2
         << setw(width) << setprecision(prec) << z1
         << setw(width) << setprecision(prec) << z2
         << setw(width) << setprecision(prec) << vx1
         << setw(width) << setprecision(prec) << vx2
         << setw(width) << setprecision(prec) << vy1
         << setw(width) << setprecision(prec) << vy2
         << setw(width) << setprecision(prec) << vz1
         << setw(width) << setprecision(prec) << vz2
         << endl;
  }

  return 0;
}

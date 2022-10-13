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

  arma::vec r_1 = arma::vec{20., 0., 20.};
  arma::vec r_2 = arma::vec{25., 25., 0.};
  arma::vec v_1 = arma::vec{0., 25., 0.};
  arma::vec v_2 = arma::vec{0., 40., 5.};

  Particle p_1 = Particle(q, m, r_1, v_1);
  Particle p_2 = Particle(q, m, r_2, v_2);

  trap.add_particle(p_1);
  trap.add_particle(p_2);

  int total_time = 50;   // [mu*s]
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

  cout << setw(width) << t(0)
       << setw(width) << r_1(0)
       << setw(width) << r_2(0)
       << setw(width) << r_1(1)
       << setw(width) << r_2(1)
       << setw(width) << r_1(2)
       << setw(width) << r_2(2)
       << setw(width) << v_1(0)
       << setw(width) << v_2(0)
       << setw(width) << v_1(1)
       << setw(width) << v_2(1)
       << setw(width) << v_1(2)
       << setw(width) << v_2(2)
       << endl;

  for (int i = 1; i < steps; i++)
  {
    trap.evolve_RK4(dt, false);

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

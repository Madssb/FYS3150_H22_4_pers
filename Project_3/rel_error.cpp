// Code for testing relative error
#include "particle.hpp"
#include "penning_trap.hpp"
#include "analytical.hpp"
#include <cmath>

using namespace std;

int main()
{
  int width = 15;
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

  double total_time = 1.;

  cout << "#" << setw(width) << "dt"
       << setw(width) << "time"
       << setw(width) << "r_ana"
       << setw(width) << "r_num"
       << setw(width) << "rel_err"
       << endl;

  for (int i = 1; i < 6; i++)
  {
    double dt = pow(10, -i);
    int steps = total_time / dt;
    arma::vec t = arma::linspace(0, total_time, steps);
    trap.particles[0].r = r;

    for (int j = 0; j < steps; j++)
    {
      trap.evolve_RK4(dt, true);

      arma::vec r_num = trap.particles[0].r;
      arma::vec r_ana = analytical_solution(r, v, p.q, B_0, p.m, V_0, d, t(j));

      double rel_err = arma::norm(arma::abs(r_num - r_ana) / arma::abs(r_ana));

      cout << setw(width) << setprecision(prec) << dt
           << setw(width) << setprecision(prec) << t(j)
           << setw(width) << setprecision(prec) << arma::norm(r_ana)
           << setw(width) << setprecision(prec) << arma::norm(r_num)
           << setw(width) << setprecision(prec) << rel_err
           << endl;
    }
  }

  return 0;
}

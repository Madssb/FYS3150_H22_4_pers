// Code for testing relative error
#include "particle.hpp"
#include "penning_trap.hpp"
#include "analytical.hpp"
#include <cmath>
#include <fstream>

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

  double total_time = 10;
  arma::vec dMax = arma::vec(5);
  arma::vec h = arma::vec(5);

  ofstream outfile;
  outfile.open("data_relerr.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width) << "dt"
          << setw(width) << "time"
          << setw(width) << "r_ana"
          << setw(width) << "r_num"
          << setw(width) << "rel_err"
          << endl;

  for (int i = 1; i < 6; i++)
  {
    double dt = pow(10, -i);
    int steps = total_time / dt;

    h(i-1) = dt;

    arma::vec t = arma::linspace(0, total_time, steps);

    trap.particles[0].r = r;

    outfile << setw(width) << setprecision(prec) << dt
            << setw(width) << setprecision(prec) << t(0)
            << setw(width) << setprecision(prec) << arma::norm(r)
            << setw(width) << setprecision(prec) << arma::norm(r)
            << setw(width) << setprecision(prec) << 0
            << endl;

    for (int j = 1; j < steps; j++)
    {
      trap.evolve_RK4(dt, true);

      arma::vec r_num = trap.particles[0].r;
      arma::vec r_ana = analytical_solution(r, v, p.q, B_0, p.m, V_0, d, t(j));

      double abs_err = arma::norm(r_num - r_ana);
      double rel_err = abs_err / arma::norm(r_ana);

      if (abs_err > dMax(i - 1))
      {
        dMax(i - 1) = abs_err;
      }

      outfile << setw(width) << setprecision(prec) << dt
              << setw(width) << setprecision(prec) << t(j)
              << setw(width) << setprecision(prec) << arma::norm(r_ana)
              << setw(width) << setprecision(prec) << arma::norm(r_num)
              << setw(width) << setprecision(prec) << rel_err
              << endl;
    }
  }
  outfile.close();

  double r_err = 0.;

  for (int k = 1; k < 5; k++)
  {
    r_err += 1. / 4 * (log(dMax(k) / dMax(k-1)) / log(h(k) / h(k-1)));
  }

  cout << "Convergence rate for RK4: " << setprecision(4) << r_err << endl;

  /* From terminal:
  > Convergence rate for FE: -0.1444
  > Convergence rate for RK4: -0.04394
  */

  return 0;
}

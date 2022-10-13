// Main code of project 3
#include "particle.hpp"
#include "penning_trap_2.hpp"
#include <cmath>
#include <fstream>

using namespace std;

int main()
{
  /// Defining constants
  const double T = 9.64852558 * 10;     // [u/(mu*s)e ]
  const double V = 9.64852558 * 1e7;    // [u(mu*m)^2/(mu*s)^2e]
  const double B_0 = 1. * T;            // Magnetic field constant
  const double V_0 = 25e-3 * V;        // Electric field constant
  const double d = 500;                 // [mu*m]
  const double mCa_plus = 40.078;       // Ca+ ion mass [u]
  const double q = 1.;                  // Charge of +1 e

  arma::vec f = arma::vec{.1, .4, .7};  // Amplitudes
  double w_V = .2;                      // Frequency [MHz]

  PenningTrap trap = PenningTrap(B_0, V_0, d, f(0), w_V);

  // Setting seed for RNG
  arma::arma_rng::set_seed_random();

  int total_time = 500;
  double dt = 1.3e-2;
  double time_steps = ceil(total_time / dt);

  // double maxFreq = 0.7;
  // double minFreq = 0.5;
  double maxFreq = 2.5;
  double minFreq = .2;
  double dw = .02;
  double freq_steps = ceil((maxFreq - minFreq) / dw);

  arma::vec t = arma::linspace(0, total_time, time_steps);
  arma::vec freq_list = arma::linspace(minFreq, maxFreq, freq_steps);
  arma::mat N_particles = arma::mat(freq_steps, 3);

  // Running simulation for each value of f and w_V
  for (int i = 0; i < f.size(); i++)
  {
    // Setting the amplitude for the trap
    trap.f = f(i);

    for (int k = 0; k < freq_steps; k++)
    {
      // Removing all particles
      trap.particles.clear();

      // Filling trap
      for (int j = 0; j < 100; j++)
      {
        arma::vec r = arma::vec(3).randn() * .1 * trap.d;
        arma::vec v = arma::vec(3).randn() * .1 * trap.d;

        Particle p = Particle(q, mCa_plus, r, v);
        trap.add_particle(p);
      }

      // Starting w/ 100 particles
      N_particles.col(i)(0) = trap.particles.size();

      // Setting the frequency
      trap.w_V = freq_list(k);

      for (int l = 0; l < time_steps; l++)
      {
        trap.evolve_RK4(t(l), dt, false);
      }

      if (k != 0)
      {
        N_particles.col(i)(k) = trap.count_particles();
      }
    }
  }

  // Getting the fraction of remaining particles
  N_particles /= 100;

  int width = 15;
  int prec = 6;

  ofstream outfile;
  outfile.open("no_of_particles.txt", ofstream::out | ofstream::trunc);

  outfile << "#" << setw(width - 1) << "freq"
                 << setw(width) << "f1"
                 << setw(width) << "f2"
                 << setw(width) << "f3"
                 << setw(width) << "N1 left"
                 << setw(width) << "N2 left"
                 << setw(width) << "N3 left"
                 << endl;

  for (int i = 0; i < freq_steps; i++)
  {
    outfile << setw(width) << setprecision(prec) << freq_list(i)
            << setw(width) << setprecision(prec) << f(0)
            << setw(width) << setprecision(prec) << f(1)
            << setw(width) << setprecision(prec) << f(2)
            << setw(width) << setprecision(prec) << N_particles.col(0)(i)
            << setw(width) << setprecision(prec) << N_particles.col(1)(i)
            << setw(width) << setprecision(prec) << N_particles.col(2)(i)
            << endl;
  }

  outfile.close();

  return 0;
}

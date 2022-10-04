// Code containing the PenningTrap class definition
#include "penning_trap.hpp"
#include "particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in)
{
  // Assigning member variables
  B_0 = B0_in;
  V_0 = V0_in;
  d = d_in;
}

// Method to add particle to trap
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}

// Method to calculate ext. E-field
arma::vec PenningTrap::external_E_field(arma::vec r)
{
  const double A = V_0 / (d * d);
  arma::vec E = arma::vec{r(0), r(1), -2. * r(2)};

  return E;
}

// Method to calculate ext. B-field
arma::vec PenningTrap::external_B_field(arma::vec r)
{
  arma::vec B = arma::vec{0., 0., B_0};

  return B;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(int i, int j)
{
  const double k_e = 1.38935333 * 1e5 ;  // Electric constant in u(mu*m)^3/(mu*s)^2e^2

  Particle p_i = particles[i];
  Particle p_j = particles[j];

  double q_j = p_j.q;
  arma::vec r_i = p_i.r;
  arma::vec r_j = p_j.r;

  arma::vec E = k_e * p_j.q * (r_i - r_j) /
                arma::pow(arma::abs(r_i - r_j), 3.0);

  return E;
}

// Total force on particle_i from external fields
arma::vec PenningTrap::total_force_external(int i)
{
  Particle p = particles[i];
  arma::vec r = p.r;

  arma::vec F = PenningTrap::external_E_field(r)
              + PenningTrap::external_B_field(r);

  return F;
}

// Total force on particle_i from other particles
arma::vec PenningTrap::total_force_particles(int i)
{
  Particle p = particles[i];
  int N = particles.size();

  arma::vec F_tot = arma::vec(3);

  for (int j = 0; j < N; j++)
  {
    if (j == i)
    {}

    else
    {
      F_tot += PenningTrap::force_particle(i, j);
    }
  }
  return F_tot;
}

// Total force on particle_i from external field and other particles
arma::vec PenningTrap::total_force(int i)
{
  arma::vec F_tot = PenningTrap::total_force_external(i)
                  + PenningTrap::total_force_particles(i);

  return F_tot;
}

// Evolve the system one time step using Runge-Kutta 4th order
// void evolve_RK4(double dt)
// {
//
// }

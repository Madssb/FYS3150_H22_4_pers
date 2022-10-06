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
arma::vec PenningTrap::external_B_field(arma::vec v)
{
  double vx = v(0), vy = v(1);

  arma::vec B = B_0 * arma::vec{vy, - vx, 0.};

  return B;
}

// Force on particle from all other particles
arma::vec PenningTrap::total_force_particles(int i, arma::vec r)
{
  // Electric constant in u(mu*m)^3/(mu*s)^2e^2
  const double k_e = 1.38935333 * 1e5;

  int N = particles.size();

  arma::vec E = arma::vec(3);

  for (int i = 0; i < N; i++)
  {
    Particle p_j = particles[i];
    double q_j = p_j.q;
    arma::vec r_j = p_j.r;

    if (i == j)
    {}

    else
    {
      E += k_e * q_j * (r - r_j) / arma::pow(arma::abs(r - r_j), 3.0);
    }
  }
  return E;
}

// Force on particle_i from particle_j
/*
arma::vec PenningTrap::force_particle(int i, int j)
{
  const double k_e = 1.38935333 * 1e5;  // Electric constant in u(mu*m)^3/(mu*s)^2e^2

  Particle p_i = particles[i];
  Particle p_j = particles[j];

  double q_j = p_j.q;
  arma::vec r_i = p_i.r;
  arma::vec r_j = p_j.r;

  arma::vec E = k_e * p_j.q * (r_i - r_j) /
                arma::pow(arma::abs(r_i - r_j), 3.0);

  return E;
}
*/

// Total force on particle
arma::vec PenningTrap::total_force(int i, double q, arma::vec r,
                                    arma::vec v)
{
  arma::vec E_ext = PenningTrap::external_E_field(r);
  arma::vec E_particles = PenningTrap(i, r)
  arma::vec E = E_ext + E_particles;

  arma::vec B = PenningTrap::external_B_field(v);

  arma::vec F_tot = q * E + q * B;

  return F_tot;
}

// Total force on particle_i from other particles
/*
arma::vec PenningTrap::total_force_particles(int i)
{
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
*/

// Total force on particle_i from external field and other particles
/*
arma::vec PenningTrap::total_force(int i)
{
  arma::vec F_tot = PenningTrap::total_force_external(i)
                  + PenningTrap::total_force_particles(i);

  return F_tot;
}
*/

// Evolve the system one time step using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double dt)
{
  int N = particles.size();
  double m;

  arma::vec a, v, r;
  arma::vec K1_r, K2_r, K3_r, K4_r;
  arma::vec K1_v, K2_v, K3_v, K4_v;

  for (int i = 0; i < N; i++)
  {
    Particle p = particles[i];

    m = p.m;
    v = p.v;
    r = p.r;
    a = PenningTrap::total_force(i) / m;

    K1_v = dt * a;
    K2_v = dt * (a + .5 * K1_v);
    K3_v = dt * (a + .5 * K2_v);
    K4_v = dt * (a + K3_v);

    v += 1. / 6 * (K1_v + 2 * K2_v + 2 * K3_v + K4_v);

    K1_r = dt * v;
    K2_r = dt * (v + .5 * K1_r);
    K3_r = dt * (v + .5 * K2_r);
    K4_r = dt * (v + K3_r);

    r += 1. / 6 * (K1_r + 2 * K2_r + 2 * K3_r + K4_r);
  }
}

// Evolve the system one time step using Forward Euler
void PenningTrap::evolve_FE(double dt)
{
  int N = particles.size();
  double m;

  arma::vec a, v, r;

  for (int i = 0; i < N; i++)
  {
    Particle p = particles[i];

    m = p.m;
    v = p.v;
    r = p.r;
    a = PenningTrap::total_force(i) / m;

    r += dt * v;
    v += dt * a;
  }
}

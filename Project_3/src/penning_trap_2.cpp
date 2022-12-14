// Code containing the PenningTrap class definition
#include "penning_trap_2.hpp"
#include "particle.hpp"

// Constructor
PenningTrap::PenningTrap(double B0_in, double V0_in, double d_in, double f_in,
                         double w_V_in)
{
  // Assigning member variables
  B_0 = B0_in;
  V_0 = V0_in;
  d = d_in;
  f = f_in;
  w_V = w_V_in;
}

// Method to add particle to trap
void PenningTrap::add_particle(Particle p_in)
{
  particles.push_back(p_in);
}

// Method to calculate ext. E-field
arma::mat PenningTrap::external_E_field(double t, arma::mat R)
{
  int N = particles.size();
  arma::mat E = arma::mat(3, N);

  // Time dependent electric field.
  double V_0_t = V_0 * (1. + f * cos(w_V * t));

  double A = V_0_t / (d * d);

  for (int i = 0; i < N; i++)
  {
    arma::vec r = R.col(i);

    // Making sure the particle escapes the trap if
    // it's distance from the center is greater than the trap.
    if (arma::norm(r) <= d)
    {
      E.col(i) = A * arma::vec{r(0), r(1), -2. * r(2)};
    }

    else
    {
      E.col(i) = arma::vec{0., 0., 0.};
    }
  }
  return E;
}

// Method to calculate ext. B-field
arma::mat PenningTrap::external_B_field(arma::mat R, arma::mat V)
{
  int N = particles.size();
  arma::vec B = arma::vec{0., 0., B_0};
  arma::mat B_tot = arma::mat(3, N);

  for (int i = 0; i < N; i++)
  {
    arma::vec r_i = R.col(i);

    if (arma::norm(r_i) < d)
    {
      arma::vec v_i = V.col(i);
      double q_i = particles[i].q;

      B_tot.col(i) = arma::cross(q_i * v_i, B);
    }

    else
    {
      B_tot.col(i) = arma::vec{0., 0., 0.};
    }
  }
  return B_tot;
}

// Force on particle_i from particle_j
arma::vec PenningTrap::force_particle(arma::vec r_i, arma::vec r_j,
                                      double q_i, double q_j)
{
  // Electric constant in u(mu*m)^3/(mu*s)^2e^2
  const double k_e = 1.38935333 * 1e5;

  arma::vec E = k_e * q_j * (r_i - r_j) / std::pow(arma::norm(r_i - r_j), 3.0);

  return q_i * E;
}

// // Total force on particle_i from other particles
arma::mat PenningTrap::total_force_particles(arma::mat R)
{
  int N = particles.size();

  arma::mat F_tot = arma::mat(3, N);
  arma::vec Fij;

  for (int i = 0; i < N; i++)
  {
    arma::vec r_i = R.col(i);
    double q_i = particles[i].q;

    // Making use of Newton's 3rd law; The force on particle i from particle j
    // is the same as the force on particle j from particle i, only with
    // opposite sign.
    for (int j = 0; j < i - 1; j++)
    {
      double q_j = particles[j].q;
      arma::vec r_j = R.col(j);

      if (arma::norm(r_i - r_j) < .3 * d)
      {
        Fij = force_particle(r_i, r_j, q_i, q_j);
      }

      else
      {
        Fij = arma::vec{0, 0, 0};
      }

      F_tot.col(i) += Fij;
      F_tot.col(j) -= Fij;
    }
  }
    return F_tot;
  }

// Total ext. force
arma::mat PenningTrap::total_force_external(double t, arma::mat R, arma::mat V)
{
  arma::mat extE = external_E_field(t, R);
  arma::mat extB = external_B_field(R, V);

  arma::mat extF_tot = extE + extB;

  return extF_tot;
}

// Total force on particle_i from external field and other particles
arma::mat PenningTrap::total_force(double t, arma::mat R, arma::mat V,
                                   bool particle_interaction)
{
  arma::mat F_tot;

  // Function to include/exclude particle interaction
  if (particle_interaction)
  {
    F_tot = total_force_external(t, R, V) + total_force_particles(R);
  }

  else
  {
    F_tot = total_force_external(t, R, V);
  }

  return F_tot;
}

// Evolve the system one time step using Runge-Kutta 4th order
void PenningTrap::evolve_RK4(double t, double dt, bool particle_interaction)
{

  int N = particles.size();

  arma::mat R = arma::zeros(3, N);
  arma::mat V = arma::zeros(3, N);
  arma::mat a = arma::zeros(3, N);
  double m = particles[0].m;

  arma::mat K1_r, K2_r, K3_r, K4_r;
  arma::mat K1_v, K2_v, K3_v, K4_v;

  for (int i = 0; i < N; i++)
  {
    R.col(i) = particles[i].r;
    V.col(i) = particles[i].v;
  }

  a = 1 / m * total_force(t, R, V, particle_interaction);

  K1_v = dt * a;
  K1_r = dt * V;

  a = 1 / m * total_force(t + .5 * dt, R + .5 * K1_r, V
                            + .5 * K1_v, particle_interaction);

  K2_v = dt * a;
  K2_r = dt * (V + 0.5 * K1_v);

  a = 1 / m * total_force(t + .5 * dt, R + .5 * K2_r, V
                            + .5 * K2_v, particle_interaction);

  K3_v = dt * a;
  K3_r = dt * (V + .5 * K2_v);

  a = 1 / m * total_force(t + dt, R + K3_r, V + K3_v, particle_interaction);

  K4_v = dt * a;
  K4_r = dt * (V + K3_v);

  R += 1. / 6 * (K1_r + 2 * K2_r + 2 * K3_r + K4_r);
  V += 1. / 6 * (K1_v + 2 * K2_v + 2 * K3_v + K4_v);

  for (int i = 0; i < N; i++)
  {
    particles[i].r = R.col(i);
    particles[i].v = V.col(i);
  }
}

// Evolve the system one time step using Forward Euler
void PenningTrap::evolve_FE(double t, double dt, bool particle_interaction)
{
  int N = particles.size();
  arma::mat R = arma::zeros(3, N);
  arma::mat V = arma::zeros(3, N);
  arma::mat a = arma::zeros(3, N);
  double m = particles[0].m;

  for (int i = 0; i < N; i++)
  {
    R.col(i) = particles[i].r;
    V.col(i) = particles[i].v;
  }

  a = 1 / m * total_force(t, R, V, particle_interaction);

  R += dt * V;
  V += dt * a;

  for (int i = 0; i < N; i++)
  {
    particles[i].r = R.col(i);
    particles[i].v = V.col(i);
  }
}

// Count no. of particles in trap
int PenningTrap::count_particles()
{
  int N = particles.size();

  // Looping over particles backwards in order to not miss any particles 
  for (int i = N - 1; i >= 0; i--)
  {
    if (arma::norm(particles[i].r) > d)
    {
      particles.erase(particles.begin() + i);
    }
  }

  return particles.size();
}

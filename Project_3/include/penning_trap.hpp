// Header file of the class PenningTrap
#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include <armadillo>
#include <string>
#include <vector>
#include "particle.hpp"

class PenningTrap
{
  public:

  // Member variables. mu is micro
  double B_0;                       // Magnetic field strength [u/mu*s*e]
  double V_0;                       // Applied potential [u*(mu*m)^2/(mu*s)^2*e]
  double d;                         // Characteristic dimension [mu*m]
  std::vector<Particle> particles;  // Vector containing Particle objects

  // Constructor
  PenningTrap(double B0_in, double V0_in, double d_in);

  // Add a particle to the trap
  void add_particle(Particle p_in);

  // External electric field at point r=(x,y,z)
  arma::mat external_E_field(arma::mat R);

  // External magnetic field at point r=(x,y,z)
  arma::mat external_B_field(arma::mat V);

  // Force on particle_i from particle_j
  arma::vec force_particle(arma::vec r_i, arma::vec r_j,
                           double q_i, double q_j);

  // Total force on particle_i from external fields
  arma::mat total_force_external(arma::mat R, arma::mat V);

  // Total force on particle_i from other particles
  arma::mat total_force_particles(arma::mat R);

  // Total force on particle_i from external field and other particles
  arma::mat total_force(arma::mat R, arma::mat V, bool particle_interaction);

  // Evolve the system one time step using Runge-Kutta 4th order
  void evolve_RK4(double dt, bool particle_interaction=true);

  // Evolve the system one time step using Forward Euler
  void evolve_FE(double dt, bool particle_interaction=true);
};


#endif

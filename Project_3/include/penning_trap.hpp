// Header file of the class PenningTrap
#ifndef __PenningTrap_hpp__
#define __PenningTrap_hpp__

#include<armadillo>
#include<string>
#include<vector>
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
  arma::vec external_E_field(arma::vec r);

  // External magnetic field at point r=(x,y,z)
  arma::vec external_B_field(arma::vec r);

  // Force on particle_i from particle_j
  arma::vec force_particle(int i, int j);

  // Total force on particle_i from external fields
  arma::vec total_force_external(int i);

  // Total force on particle_i from other particles
  arma::vec total_force_particles(int i);

  // Total force on particle_i from external field and other particles
  arma::vec total_force(int i);

  // Evolve the system one time step using Runge-Kutta 4th order
  void evolve_RK4(double dt);

  // Evolve the system one time step using Forward Euler
  void evolve_FE(double dt);
};


#endif

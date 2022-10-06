// Code containing the definition of class Particle
#include "particle.hpp"

// Constructor
Particle::Particle(double charge_in, double mass_in, arma::vec position_in,
arma::vec velocity_in)
{
  // Assigning member variables
  q = charge_in;
  m = mass_in;
  r = position_in;
  v = velocity_in;
}

// Method that returns particle properties
std::string Particle::properties()
{
  std::string info_string = "Particle properties:\nCharge: " + std::to_string(q)
                          + " e.\nMass: " + std::to_string(m) + " u";

  return info_string;
}

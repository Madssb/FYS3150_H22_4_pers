// Header file of the class Particle
#ifndef __Particle_hpp__
#define __Particle_hpp__

#include<armadillo>
#include<string>

class Particle
{
  public:

    // Member variables
    double q;       // Charge [e]
    double m;       // Mass [u]
    arma::vec r;    // Position [micrometer]
    arma:: vec v;   // Velocity [micrometer/microsecond]

    // Constructor
    Particle(double charge_in, double mass_in, arma::vec position_in,
    arma::vec velocity_in);

    // Method that returns particle properties
    std::string properties();
};

#endif

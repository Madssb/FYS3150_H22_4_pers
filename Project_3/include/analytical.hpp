// Contains the analytical solution i z-direction
#ifndef __Analytical_hpp__
#define __Analytical_hpp__

#include <armadillo>

/// Defining constants
const double T = 9.64852558 * 10;   // [u/(mu*s)e ]
const double V = 9.64852558 * 1e7;  // [u(mu*m)^2/(mu*s)^2e]
const double B_0 = 1. * T;          // Magnetic field constant
const double V_0 = 10. * V;         // Electric field constant
const double d = 1e4;               // [mu*m]

// Analytical solution
arma::vec analytical_solution(arma::vec r_in, arma::vec v_in, double q,
                              double B_ext, double m, double V_ext,
                              double d, double t);

#endif

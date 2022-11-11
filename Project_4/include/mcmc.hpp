/*
Header file for the Markov Chain Monte Carlo (MCMC) method
*/

#ifndef __MCMC_hpp__
#define __MCMC_hpp__

#include <armadillo>

// Perform one MCMC cycle, ie. N attempted spin flips
void mcmc(arma::mat& lattice, std::mt19937& generator, int L, double& E,
          double& M, double T);

#endif

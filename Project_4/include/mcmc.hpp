/*
This is the header file of the Markov Chain Monte Carlo (MCMC) method
*/

#ifndef __MCMC_hpp__
#define __MCMC_hpp__

#include <armadillo>

// Doing a full MCMC cycle to the lattice
void mcmc(arma::mat& lattice, int L, double T, std::mt19937& generator);

#endif

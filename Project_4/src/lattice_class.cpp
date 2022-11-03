/*
Definition of the class Lattice
*/

#include "lattice_class.hpp"

// Constructor
Lattice::Lattice(int L_in, double T_in, int nCycles_in)
{
  // Assigning member variables
  L = L_in;
  T = T_in;
  nCycles = nCycles_in;
}

// Creating lattice
arma::mat Lattice::create_spins()
{
  // Setting seed of RNG
  arma::arma_rng::set_seed_random();

  // Filling lattice with random spins
  arma::mat lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0, 1));

  // Setting all 0s to -1s
  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      if (lattice(i, j) == 0)
      {
        lattice(i, j) = -1;
      }
    }
  }

  return lattice;
}

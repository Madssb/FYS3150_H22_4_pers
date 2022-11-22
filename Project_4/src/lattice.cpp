/*
Initializing the lattice
*/

#include "lattice.hpp"

arma::mat initialize_lattice(int L, double temp, double& E, double& M,
                             bool ordered)
{
  arma::mat lattice;

  if (ordered)
  {
    lattice = arma::mat(L, L, arma::fill::ones);
  }

  else
  {
    // Setting RNG seed
    arma::arma_rng::set_seed_random();

    // Filling lattice with random spins
    lattice = arma::randi<arma::mat>(L, L, arma::distr_param(0, 1));

    // Setting all 0s to -1s
    for (int x = 0; x < L; x++)
    {
      for (int y = 0; y < L; y++)
      {
        if (lattice(x, y) == 0)
        {
          lattice(x, y) = -1;
        }
      }
    }
  }

  // Computing initial energy and magnetization
  for (int x = 0; x < L; x++)
  {
    for (int y = 0; y < L; y++)
    {
      double spin = lattice(x, y);
      double up = lattice((x + 1) % L, y);
      double down = lattice(x, (y + 1) % L);

      E -= spin * (up + down);
      M += spin;
    }
  }

  return lattice;
}

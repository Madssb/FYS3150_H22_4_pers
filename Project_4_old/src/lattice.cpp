/*
Definition of the initialization and computing functions of the lattice
*/

#include "lattice.hpp"

// Creating lattice
arma::mat create_lattice(int L)
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

// Computing energy of spin (i, j)
int energySpin_ij(arma::mat& lattice, int L, int i, int j)
{
  int spin = lattice(i, j);
  int spinLeft = lattice(i, (j - 1 + L) % L);
  int spinRight = lattice(i, (j + 1) % L);
  int spinUp = lattice((i - 1 + L) % L, j);
  int spinDown = lattice((i + 1) % L, j);

  int E = -spin * (spinLeft + spinRight + spinUp + spinDown);

  return E;
}

// Computing magnetization of lattice
int magnetizationLattice(arma::mat& lattice, int L)
{
  int M = 0;

  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      M += lattice(i, j);
    }
  }

  return M;
}

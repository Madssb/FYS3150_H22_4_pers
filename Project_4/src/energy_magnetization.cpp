/*
Definition of the energy and magnetization functions
*/

#include "energy_magnetization.hpp"

int totE(arma::mat& lattice, int L)
{
  int E = 0;

  // Given the periodic boundary condition, we can multiply spin (i, j) with
  // the spin to the right and below. In this way we will avoid double counting.
  for (int i = 0; i < L; i++)
  {
    for (int j = 0; j < L; j++)
    {
      int spin = lattice(i, j);
      int spinRight = lattice(i, (j + 1) % L);
      int spinDown = lattice((i + 1) % L, j);
      E -= spin * (spinRight + spinDown);
    }
  }

  return E;
}

int magnetizationLattice(arma::mat& lattice)
{
  int M = arma::accu(lattice);

  return M;
}

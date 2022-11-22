/*
Perform one MCMC cycle
*/

#include "mcmc.hpp"
#include <random>
#include <map>
#include <cmath>

void mcmc(arma::mat& lattice, std::mt19937& generator, int L, double& E,
          double& M, double T)
{
  // Construct a distribution of uniformly random integers in [0, L-1]
  std::uniform_int_distribution<int> randInt(0, L - 1);

  // Same as above for floats in [0.0, 1.0)
  std::uniform_real_distribution<double> randFloat(0., 1.);

  // For "nsteps" number of times, we select a spin, flip it, calculating the
  // energy difference (dE) and accepting/rejecting. We know that there are
  // only five possible values dE can take, so we will use that.
  std::map<int, double> dE;
  dE[8] = std::exp(-8 / T);
  dE[4] = std::exp(-4 / T);
  dE[0] = 1.;
  dE[-4] = 1.;
  dE[-8] = 1.;

  for (int x = 0; x < L; x++)
  {
    for(int y = 0; y < L; y++)
    {
      int ix = randInt(generator);
      int iy = randInt(generator);

      int spin = lattice(ix, iy);
      int up = lattice((ix + L - 1) % L, iy);
      int down = lattice((ix + 1) % L, iy);
      int left = lattice(ix, (iy + L - 1) % L);
      int right = lattice(ix, (iy + 1) % L);

      int deltaE = 2 * spin * (up + down + left + right);
      double r = randFloat(generator);

      double A = dE[deltaE];

      if(r < A)
      {
        lattice(ix, iy) *= -1;

        M += (double) 2 * lattice(ix, iy);
        E += (double) deltaE;
      }
    }
  }
}

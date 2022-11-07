/*
Definition of the MCMC method.
*/

#include "mcmc.hpp"
#include "lattice.hpp"
#include <random>
#include <chrono>
#include <map>
#include <cmath>

// Doing a full MCMC cycle by randomly selecting a spin, flipping it and
// accept/reject according to the Metropolis algorithm
int mcmc(arma::mat& lattice, int L, double T, std::mt19937& generator)
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

  int N = L * L;
  int diffE = 0;

  // Making a copy of the lattice and flipping the spin
  arma::mat testLattice = lattice;

  // Doing a full MC cycle
  for (int n = 0; n < N; n++)
  {
    int i = randInt(generator);
    int j = randInt(generator);
    int spin = lattice(i, j);
    int E_0 = energySpin_ij(lattice, L, i, j);

    testLattice(i, j) = -spin;

    int E_1 = energySpin_ij(testLattice, L, i, j);
    int energyDiff = E_1 - E_0;

    // Checking if energyDiff is accepted or rejected
    double r = randFloat(generator);
    int A = dE[energyDiff];

    if (r < A)
    {
      lattice = testLattice;
      diffE += energyDiff;
    }

    else
    {
      testLattice(i, j) = spin;
    }
  }

  return diffE;
}

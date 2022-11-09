/*
Source file for the functions to calculate expected values
*/

#include "expectation_values.hpp"
#include <map>
#include <iostream>
#include <cmath>

// Calculate the expected energy
std::vector<double> expected_energy(std::vector<double> epsilonEnergies,
                                    double initEnergy, double T)
{
  double N = epsilonEnergies.size();
  std::vector<double> expEpsilon;
  double expectedEnergy = initEnergy;

  // Partition funct. for 2x2 case
  double Z = 4 * (std::cosh(8 / T) + 3);

  // Probabilities of energies in the 2x2 case
  std::map<int, double> probs;
  std::cout << expectedEnergy << std::endl;
  probs[-2] = 1. / Z * std::exp(8 / T);
  probs[2] = 1. / Z * std::exp(-8 / T);
  probs[0] = 1. / Z;

  for (int n = 0; n < N; n++)
  {
    int E = epsilonEnergies[n];
    expectedEnergy += (E * probs[E]);
    expEpsilon.push_back(expectedEnergy / (N - 1));
  }

  return expEpsilon;
}
